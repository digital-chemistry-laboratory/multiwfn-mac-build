!======================================================================================
!======================================================================================
!---------- Study various real space functions at a point (i.e. 0 dimension) ----------
!======================================================================================
!======================================================================================
subroutine study0dim
use defvar
implicit real*8 (a-h,o-z)
character inpstring*80
integer :: iprintfunc=1 !The default function whose gradient and Hessian will be outputted at a point by main function 1
do while(.true.)
	write(*,*)
	write(*,*) "          ------------ Calculate properties at a point ------------ "
	write(*,"(a)") " Now input X,Y,Z of the point to be studied in Bohr or Angstrom, e.g. 3.3,2.0,-0.3"
	write(*,"(a)") " or input e.g. ""a5"" to use nuclear position of atom 5"
	write(*,"(a)") "    input e.g. ""o8"" to select orbital 8, whose wavefunction value will be shown"
	write(*,"(a)") "    input e.g. ""f3"" to select function 3, whose gradient and Hessian will be shown, input ""allf"" can print all available functions"	
	write(*,"(a)") "    input ""d"" can decompose a real space function at specified point to orbital contributions"
	write(*,"(a)") "    input ""q"" can return"
	read(*,"(a)") inpstring
	inpstring=adjustl(inpstring)
	if (inpstring(1:1)=='q') then
		return
	else if (inpstring(1:4)=='allf') then
		call funclist
	else if (inpstring(1:1)=='f') then	
		read(inpstring(2:),*) iprintfunc
		write(*,"(' OK, real space function',i5,' has been selected')") iprintfunc
	else if (inpstring(1:1)=='a') then
		read(inpstring(2:),*) iatm
		if (iatm>0.and.iatm<=ncenter) then
			iskipnuc_old=iskipnuc
			iskipnuc=iatm
			write(*,*) "Note: Unless otherwise specified, all units are in a.u."
			write(*,"(' Atom',i6,'  X,Y,Z(Bohr):',3f14.8)") iatm,a(iatm)%x,a(iatm)%y,a(iatm)%z
			call showptprop(a(iatm)%x,a(iatm)%y,a(iatm)%z,iprintfunc,6)
			iskipnuc=iskipnuc_old
		else
			write(*,"(' The atomic index exceeded valid range, should be <=',i5)") ncenter
		end if
	else if (inpstring(1:1)=='o') then
		read(inpstring(2:),*) iorbseltmp
		if (iorbseltmp>0.and.iorbseltmp<=nmo) then
			iorbsel=iorbseltmp
		else
			write(*,"(' The orbital index exceeded valid range, should be larger than 0 and <=',i5)") nmo
		end if
		write(*,"(' OK, orbital',i5,' has been selected')") iorbsel
	else if (inpstring(1:1)=='d') then
		write(*,*) "Input X,Y,Z of the point, e.g. 3.3,2.0,-0.3"
		read(*,*) xin,yin,zin
		write(*,*) "You inputted coordinate is in which unit?  1: Bohr  2: Angstrom"
		read(*,*) iunit
		if (iunit==2) then
			xin=xin/b2a
			yin=yin/b2a
			zin=zin/b2a
		end if
		call decompptprop(xin,yin,zin)
	else
		read(inpstring,*,iostat=ierror) xin,yin,zin
        if (ierror/=0) then
			write(*,*) "Error: Unable to recognize your input!"
			cycle
        end if
		write(*,*) "You inputted coordinate is in which unit?  1: Bohr  2: Angstrom"
		read(*,*) iunit
		if (iunit==2) then
			xin=xin/b2a
			yin=yin/b2a
			zin=zin/b2a
		end if
		write(*,*) "Note: Unless otherwise specified, all units are in a.u."
		call showptprop(xin,yin,zin,iprintfunc,6)
	end if
end do
end subroutine










!======================================================================================
!======================================================================================
!-------- Study various real space functions along a line (i.e. 1 dimension) ----------
!======================================================================================
!======================================================================================
subroutine study1dim
use defvar
use functions
use plot
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200,c400tmp*400,filename_tmp*200,graphformat_old*4

!Clean custom operation setting (global)that possibly defined by other modules
ncustommap=0
if (allocated(custommapname)) deallocate(custommapname)
if (allocated(customop)) deallocate(customop)
ipromol=0

do while(.true.)
	write(*,*) "-10 Return to main menu"
	if (ncustommap==0) then
		write(*,*) "-2 Obtain deformation property"
		write(*,*) "-1 Obtain promolecule property"
		write(*,*) "0 Set custom operation"
	end if
	call selfunc_interface(1,isel) !Interface of selecting real space function
	
	if (isel==-10) then
		return
	else if (isel==-2) then
		call setPromol
		customop='-'
	else if (isel==-1) then
		ipromol=1 !Special case, obtain promolecular property
		call setPromol
		customop='+'
	else if (isel==0) then
		call customplotsetup
	else if (isel==111) then
		write(*,*) "Input indices of two atoms to calculate Becke overlap weight, e.g. 1,4"
		write(*,*) "or input index of an atom and zero to calculate Becke atomic weight, e.g. 5,0"
		read(*,*) iatmbecke1,iatmbecke2
	else
		exit
	end if
end do

do while(.true.)
	write(*,"(a,f8.4,a)") " 0 Set extension distance for mode 1, current:",aug1D," Bohr"
	write(*,*) "1 Input index of two atoms to define a line"
	write(*,*) "2 Input coordinate of two points to define a line"
	read(*,*) isel2

	if (isel2==0) then
		write(*,*) "Input extension distance (in Bohr, e.g. 2.5)"
		read(*,*) aug1D
	else if (isel2==1) then
		do while(.true.)
			write(*,*) "Input two indices to select two atoms, e.g. 1,3"
			read(*,*) iselatm1,iselatm2
			if (iselatm1/=iselatm2.and.min(iselatm1,iselatm2)>=1.and.max(iselatm1,iselatm2)<=ncenter) exit
			write(*,*) "Error: Invalid input!"
		end do
		write(*,"(' Atom',i5,'  Charge:',f6.2,'  X,Y,Z:',3f12.6,' Bohr')") iselatm1,a(iselatm1)%charge,a(iselatm1)%x,a(iselatm1)%y,a(iselatm1)%z
		write(*,"(' Atom',i5,'  Charge:',f6.2,'  X,Y,Z:',3f12.6,' Bohr')") iselatm2,a(iselatm2)%charge,a(iselatm2)%x,a(iselatm2)%y,a(iselatm2)%z
		torgx=a(iselatm1)%x
		torgy=a(iselatm1)%y
		torgz=a(iselatm1)%z
		tendx=a(iselatm2)%x
		tendy=a(iselatm2)%y
		tendz=a(iselatm2)%z
		ratio=dsqrt((tendx-torgx)**2+(tendy-torgy)**2+(tendz-torgz)**2)/aug1D
		orgx1D=torgx-(tendx-torgx)/ratio
		orgy1D=torgy-(tendy-torgy)/ratio
		orgz1D=torgz-(tendz-torgz)/ratio
		endx1D=tendx+(tendx-torgx)/ratio
		endy1D=tendy+(tendy-torgy)/ratio
		endz1D=tendz+(tendz-torgz)/ratio
		totdist=dsqrt((endx1D-orgx1D)**2+(endy1D-orgy1D)**2+(endz1D-orgz1D)**2)
		atomr1=aug1D
		atomr2=totdist-aug1D
		exit
	else if (isel2==2) then
		write(*,*) "Input x1,y1,z1,x2,y2,z2 to define two points (in Bohr)"
		write(*,*) "e.g. 0.5,0,3.2,-1,-0.26,2.8"
        write(*,*) "To input coorindates in Angstrom, you can add "" A"" suffix, e.g. 0,0,0,0,0,5 A"
        read(*,"(a)") c80tmp
		read(c80tmp,*) x1,y1,z1,x2,y2,z2
        if (index(c80tmp,'A')/=0) then
            x1=x1/b2a;y1=y1/b2a;z1=z1/b2a;x2=x2/b2a;y2=y2/b2a;z2=z2/b2a
        end if
		orgx1D=x1
		orgy1D=y1
		orgz1D=z1
		endx1D=x2
		endy1D=y2
		endz1D=z2
		atomr1=0D0
		atomr2=0D0
		exit
	end if
end do

if (isel==12) then !Calculate ESP is time consuming, so decrease the number of points than default
	npointcurve=num1Dpoints/6
	write(*,*) "Please wait..."
else
	npointcurve=num1Dpoints !The number of data to plot. num1Dpoints is defined in settings.ini
    if (ifPBC>0) then
	    if (expcutoff_PBC<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' will be ignored')") expcutoff_PBC
    else
    	if (expcutoff<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' will be ignored')") expcutoff
    end if
end if

if (allocated(curvex)) deallocate(curvex)
if (allocated(curvey)) deallocate(curvey)
if (allocated(curveytmp)) deallocate(curveytmp)
allocate(curvex(npointcurve))
allocate(curvey(npointcurve))
if (ncustommap/=0) allocate(curveytmp(npointcurve))
transx=(endx1D-orgx1D)/npointcurve
transy=(endy1D-orgy1D)/npointcurve
transz=(endz1D-orgz1D)/npointcurve
transr=dsqrt(transx**2+transy**2+transz**2)

write(*,*)
write(*,"(' Original point in X,Y,Z:    ',3f10.5,' Bohr')") orgx1D,orgy1D,orgz1D !Output grid information
write(*,"(' End point in X,Y,Z:         ',3f10.5,' Bohr')") endx1D,endy1D,endz1D
write(*,"(' Translation vector in X,Y,Z:',3f10.5,' Norm:',f9.5,' Bohr')") transx,transy,transz,transr
write(*,"(' Number of points:',i10)") npointcurve
icustom=0

if (isel==100.and.(iuserfunc==57.or.iuserfunc==58.or.iuserfunc==59)) then !Calculate g1,g2,g3 terms defined by Shubin, they rely on rho_0
    call g1g2g3line(orgx1D,orgy1D,orgz1D,transx,transy,transz)
else !Common case
    do while(.true.)
	    if (.not.(ipromol==1.and.icustom==0)) then !Calculate property for whole system. However for promolecular case at initial step, skip it
		    
            !Try to invoke cubegen to calculate ESP
		    alive=.false.
		    if (cubegenpath/=" ".and.ifiletype==1.and.isel==12) then
			    inquire(file=cubegenpath,exist=alive)
			    if (.not.alive.and.cubegenpath/="none") then
				    write(*,"(a)") " Note: Albeit current file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been defined, &
				    &the cubegen cannot be found, therefore electrostatic potential will still be calculated using internal code of Multiwfn"
			    end if
		    end if
		    if (alive.and.ifiletype==1.and.isel==12) then !Indeed use cubegen to calculate ESP
			    write(*,"(a)") " Since the input file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been properly defined, &
			    &now Multiwfn directly invokes cubegen to calculate electrostatic potential"
			
			    !Generate cubegen input file
			    open(10,file="cubegenpt.txt",status="replace")
			    do ipt=1,npointcurve
				    rnowx=orgx1D+(ipt-1)*transx
				    rnowy=orgy1D+(ipt-1)*transy
				    rnowz=orgz1D+(ipt-1)*transz
				    curvex(ipt)=ipt*transr
				    write(10,"(3f16.8)") rnowx*b2a,rnowy*b2a,rnowz*b2a
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
			    do ipt=1,npointcurve
				    read(10,*) rnouse,rnouse,rnouse,curvey(ipt)
			    end do
			    close(10)
			
			    !Delete intermediate files
				call delfile("cubegenpt.txt ESPresult.cub nouseout")
		
		    else !Normal case, use internal code to calculate data
                nthreads_old=nthreads
                !Use fast ESP evaluation function for almost all functions that rely on ESP
                if (ifdoESP(isel).and.(iESPcode==2.or.iESPcode==3)) then
                    call doinitlibreta(1)
                    if (isys==1.and.nthreads>12) nthreads=12
                end if
                ifinish=0
			    !$OMP parallel do shared(curvex,curvey,ifinish) private(ipt,rnowx,rnowy,rnowz) num_threads(nthreads)
			    do ipt=1,npointcurve  !Calculate data for line plot
				    rnowx=orgx1D+(ipt-1)*transx
				    rnowy=orgy1D+(ipt-1)*transy
				    rnowz=orgz1D+(ipt-1)*transz
				    curvex(ipt)=ipt*transr
				    if (isel==111) then
					    curvey(ipt)=beckewei(rnowx,rnowy,rnowz,iatmbecke1,iatmbecke2)
				    else
					    curvey(ipt)=calcfuncall(isel,rnowx,rnowy,rnowz)
				    end if
                    !$OMP critical
                    ifinish=ifinish+1
                    call showprog(ifinish,npointcurve)
                    !$OMP end critical
			    end do
			    !$OMP end parallel do
                nthreads=nthreads_old
		    end if
	    end if
	
	    if (ncustommap==0) then !Normal case
		    exit
	    else !Custom operation or deformation/promolecular property
		    if (icustom==0) then !First time
			    if (ipromol==0) then !Backup the whole property as curveytmp
				    curveytmp=curvey
			    else !Initial data for promolecular property is clearly zero
				    curveytmp=0
			    end if
		    else !Not first time
			    if (customop(icustom)=='+') curveytmp=curveytmp+curvey
			    if (customop(icustom)=='-') curveytmp=curveytmp-curvey
			    if (customop(icustom)=='x'.or.customop(icustom)=='*') curveytmp=curveytmp*curvey
			    if (customop(icustom)=='/') curveytmp=curveytmp/curvey
		    end if
		    if (icustom/=ncustommap) then !Not the last time
			    icustom=icustom+1
			    filename=custommapname(icustom)
                call savePBCinfo
			    call dealloall(0)
			    write(*,"(' Loading: ',a)") trim(filename)
			    call readinfile(filename,1)
                call loadPBCinfo
			    !Input the MO index for current file. Since the MO index may be not the same as the first loaded one
			    if (isel==4) then
				    write(*,"(' Input index of the orbital to be calculated for ',a,', e.g. 3')") trim(filename)
				    read(*,*) iorbsel
			    end if
		    else !Last time
			    curvey=curveytmp
			    call dealloall(0)
			    write(*,"(' Reloading: ',a)") trim(firstfilename)
			    call readinfile(firstfilename,1)
			    exit
		    end if
	    end if
    end do
end if
    
write(*,"(/,' Minimal/Maximum value:',2E16.8)") minval(curvey),maxval(curvey)
write(*,"(' Summing up all values:',E18.8,'  Integration value:',E18.8)") sum(curvey),sum(curvey)*transr
exty=(maxval(curvey)-minval(curvey))/10
if (exty<maxval(abs(curvey))/100) exty=maxval(abs(curvey))/10 !Sometimes the difference between maximal and minimal values are too small, so do special treatment
curveymin=minval(curvey)-exty
curveymax=maxval(curvey)+exty
steplabx=maxval(curvex)/10
steplaby=exty
graphformat_old=graphformat
isel=-1
curveYname=" "

do while(.true.)
	if (isilent==0.and.isel==-1) then
		if (atomr1==atomr2) then !The line was defined by two terminal, so do not plot circles to highlight atoms
			call drawcurve(curvex,curvey,npointcurve,0D0,maxval(curvex),steplabx,curveymin,curveymax,steplaby,"show",axisnamey=curveYname)
		else !The line was defined by two atoms, plot circles to highlight atoms
			call drawcurve(curvex,curvey,npointcurve,0D0,maxval(curvex),steplabx,curveymin,curveymax,steplaby,"show",atomr1,atomr2,curveYname)
		end if
	end if
    write(*,*)
    call menutitle("Post-processing menu",10,1)
	write(*,*) "-2 Set format of exported image file, current: "//graphformat
	write(*,*) "-1 Show the graph again"
	write(*,*) "0 Return to main menu"
	write(*,*) "1 Save the graph to a file in current folder"
	write(*,*) "2 Export the data to line.txt in current folder"
	write(*,"(a,1PE14.6,a,1PE14.6)") " 3 Change range of Y axis, current: from",curveymin," to",curveymax
	if (icurve_vertlinex==0) write(*,*) "4 Draw a vertical line at specific X"
	if (icurve_vertlinex==1) write(*,*) "4 Delete the vertical line"
	write(*,"(' 5 Change the ratio of X and Y length, current:',f10.5)") curvexyratio
	write(*,"(' 6 Find the positions of local minimum and maximum')")
	write(*,"(' 7 Find the positions where function value equals to specified value')")
	if (ilog10y==0) write(*,*) "8 Use logarithmic scaling of Y axis"
	if (ilog10y==1) write(*,*) "8 Use linear scaling of Y axis"
	write(*,"(a,a)") " 9 Change the line color, current: ",trim(colorname(iclrcurve))
	if (ilog10y==0) write(*,"(a,f8.3,1PE14.5)") " 10 Set label intervals in X and Y axes, current:",steplabx,steplaby
	if (ilog10y==1) write(*,"(a,f8.3)") " 10 Set label interval in X axis, current:",steplabx
	if (ilenunit1D==1) write(*,*) "11 Change length unit of the graph from Bohr to Angstrom"
	if (ilenunit1D==2) write(*,*) "11 Change length unit of the graph from Angstrom to Bohr"
	write(*,"(a,i3)") " 12 Set width of curve line, current:",icurvethick
    write(*,"(a,2i3)") " 13 Set number of digits after decimal point of the labels on X and Y axes, current:",numdiglinex,numdigliney
    write(*,"(a,i3)") " 14 Set text size of tick labels, current:",curve_axistextsize
    write(*,"(a,i3)") " 15 Set text size of axis names, current:",curve_axisnamesize
    write(*,"(a)") " 16 Set name of Y axis"

	read(*,*) isel
    if (isel==-2) then
        call setgraphformat
	else if (isel==-1) then
		cycle
	else if (isel==0) then
        graphformat=graphformat_old
		exit
	else if (isel==1) then
		if (atomr1==atomr2) then
			call drawcurve(curvex,curvey,npointcurve,0D0,maxval(curvex),steplabx,curveymin,curveymax,steplaby,"save",axisnamey=curveYname)
		else
			call drawcurve(curvex,curvey,npointcurve,0D0,maxval(curvex),steplabx,curveymin,curveymax,steplaby,"save",atomr1,atomr2,curveYname)
		end if
		write(*,"(a,a,a)") " The graph have been saved to ",trim(graphformat)," format file with ""dislin"" prefix in current directory"
	else if (isel==2) then
        c200tmp="line.txt"
        if (iaddprefix==1) call addprefix(c200tmp)
		open(10,file=c200tmp,status="replace")
		do ipt=1,npointcurve  !Output result to file
			rnowx=orgx1D+(ipt-1)*transx
			rnowy=orgy1D+(ipt-1)*transy
			rnowz=orgz1D+(ipt-1)*transz
			write(10,"(4f12.6,1PE20.10E3)") rnowx*b2a,rnowy*b2a,rnowz*b2a,curvex(ipt)*b2a,curvey(ipt)
		end do
		close(10)
		write(*,"(a)") " Data have been exported to "//trim(c200tmp)//" in current folder"
		write(*,"(a)") " Unit is Angstrom. The first three columns are actual coordinates, the fourth column &
		&is X position in the curve graph, the fifth column is function value"
	else if (isel==3) then
		if (ilog10y==0) then
			write(*,*) "Input minimum and maximum value of Y axis, e.g. -0.1,2"
		else if (ilog10y==1) then
			write(*,*) "Input minimum and maximum value of Y axis, e.g. -1,5 means from 10^-1 to 10^5"
		end if
		read(*,*) curveymin,curveymax
	else if (isel==4) then
		if (icurve_vertlinex==0) then
			write(*,*) "Input X position in Bohr, e.g. 1.78"
			read(*,*) curve_vertlinex
			icurve_vertlinex=1
		else
			icurve_vertlinex=0
		end if
	else if (isel==5) then
		write(*,*) "Input the ratio, e.g. 0.6"
		read(*,*) curvexyratio
	else if (isel==6) then
		call showcurveminmax(npointcurve,curvex,curvey,ilenunit1D)
	else if (isel==7) then
		write(*,*) "Input a value, e.g. 0.4"
		read(*,*) specvalue
		numfind=0
		do ipt=1,npointcurve-1
			if ( (specvalue>curvey(ipt).and.specvalue<curvey(ipt+1)) .or.&
			 (specvalue<curvey(ipt).and.specvalue>curvey(ipt+1)) ) then
				!Use linear interpolation to evaluate the X position
				tmpratio=abs(specvalue-curvey(ipt))/abs(curvey(ipt+1)-curvey(ipt))
				specvaluex=curvex(ipt)+tmpratio*transr
				numfind=numfind+1
				write(*,"(' #',i4,'    X position:',f12.6,' Bohr')") numfind,specvaluex
			end if
		end do
		if (numfind==0) write(*,*) "Found nothing"
	else if (isel==8) then
		if (ilog10y==1) then
			ilog10y=0
			!Recover default limit
			curveymin=minval(curvey)-(maxval(curvey)-minval(curvey))/10
			curveymax=maxval(curvey)+(maxval(curvey)-minval(curvey))/10
		else if (ilog10y==0) then
			ilog10y=1
			write(*,*) "Input minimum and maximum value of Y axis, e.g. -1,5 means from 10^-1 to 10^5"
			read(*,*) curveymin,curveymax
		end if
	else if (isel==9) then
		write(*,*) "Use which color?"
		call selcolor(iclrcurve)
	else if (isel==10) then
		if (ilog10y==0) then 
			write(*,*) "Input the step size between the labels in X and Y axes, respectively"
			write(*,*) "e.g. 1.5,20"
			read(*,*) steplabx,steplaby
		else if (ilog10y==1) then
			write(*,*) "Input the step size between the labels in X axis, e.g. 1.5"
			read(*,*) steplabx
		end if
	else if (isel==11) then
		if (ilenunit1D==1) then
			ilenunit1D=2 !Angstrom
		else if (ilenunit1D==2) then
			ilenunit1D=1 !Bohr
		end if
	else if (isel==12) then
		write(*,*) "Input line width, e.g. 5"
		read(*,*) icurvethick
	else if (isel==13) then
		write(*,*) "Input number of digits after decimal point of the labels on X and Y axes, e.g. 1,2"
        read(*,*) numdiglinex,numdigliney
    else if (isel==14) then
        write(*,*) "Input text size of ticks on the axes, e.g. 45"
        write(*,"(' Current value:',i3)") curve_axistextsize
        read(*,*) curve_axistextsize
    else if (isel==15) then
        write(*,*) "Input text size of axis names, e.g. 40"
        write(*,"(' Current value:',i3)") curve_axisnamesize
        read(*,*) curve_axisnamesize
    else if (isel==16) then
        write(*,*) "Input text of name of Y-axis, e.g. Function value"
        read(*,"(a)") curveYname
	end if
end do

end subroutine







!======================================================================================
!======================================================================================
!---------- Study various real space functions on a plane (i.e. 2 dimension) ----------
!======================================================================================
!======================================================================================
!itask=0: Invoked by main function 4
!itask=1: Invoked by NICS-2D scan
!itask=2: Invoked by (hyper)polarizability density (hyper_polar_dens)
!itype: Currently only used by itask=2, choose type of (hyper)polarizability density (>0) or spatial contribution to (hyper)polarizability (<0)
!itype2: Currently only used by itask=2, choose direction of spatial contribution to (hyper)polarizability
subroutine study2dim(itask,itype,itype2)
use defvar
use util
use topo
use functions
use GUI
implicit real*8 (a-h,o-z)
integer,allocatable :: tmparrint(:)
integer itask,itype,itype2
character c80tmp*80,c200tmp*200,c2000tmp*2000,selectyn
real*8 vec1old(3),vec2old(3),vec1(3),vec2(3),tmpmat(3,3),tmpvec(3)
real*8,allocatable :: d1add(:,:),d1min(:,:),d2add(:,:),d2min(:,:),d1addtmp(:,:),d1mintmp(:,:),d2addtmp(:,:),d2mintmp(:,:) !Store temporary data for drawing gradient map
real*8,allocatable :: planemat_cust(:,:) !For storing temporary data of doing custom map
real*8,allocatable :: planemat_bk(:,:) !Used to backup plane data

if (itask/=2) then
	ncustommap=0 !Clean custom operation setting that possibly defined by other modules
	if (allocated(custommapname)) deallocate(custommapname)
	if (allocated(customop)) deallocate(customop)
end if
if (allocated(tmparrint)) deallocate(tmparrint)
ipromol=0

if (itask==0) then
	do while(.true.)
		write(*,*) "-10 Return to main menu"
		if (ncustommap==0) then
			write(*,*) "-2 Obtain deformation property"
			write(*,*) "-1 Obtain promolecule property"
			write(*,*) "0 Set custom operation"
		end if
		call selfunc_interface(2,ifuncsel) !Interface of selecting real space function
		if (ifuncsel==-10) then
			return
		else if (ifuncsel==-2) then
			call setPromol
			customop='-'
		else if (ifuncsel==-1) then
			ipromol=1 !Special case, obtain promolecular property
			call setPromol
			customop='+'
		else if (ifuncsel==0) then
			call customplotsetup
		else if (ifuncsel==111) then !Calculate Becke weighting function
			write(*,*) "Input indices of two atoms to calculate Becke overlap weight, e.g. 1,4"
			write(*,*) "or input index of an atom and zero to calculate Becke atomic weight, e.g. 5,0"
			read(*,*) iatmbecke1,iatmbecke2
			exit
		else if (ifuncsel==112) then !Calculate Hirshfeld weighting function
			write(*,*) "Input index of the atoms that you want to calculate Hirshfeld weight"
			write(*,*) "e.g. 2,3,7-10"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,ntmp)
			allocate(tmparrint(ntmp))
			call str2arr(c2000tmp,ntmp,tmparrint)
			write(*,"(a)") " How to generate the atomic densities that used in the calculation of Hirshfeld weight?"
			write(*,*) "1 Based on atomic .wfn files"
			write(*,*) "2 Based on built-in atomic densities (see Appendix 3 of the manual for detail)"
			read(*,*) iHirshdenstype
			if (iHirshdenstype==1) call setpromol
			exit
		else if (ifuncsel==500.or.ifuncsel==510.or.ifuncsel==511.or.ifuncsel==512) then !Calculate rho(A)*ln[rho(A)/rho0(A)], or rho(A), or rho0(A)
			call setpromol
			allocate(tmparrint(1))
			write(*,*) "Input index of the atom you are interested in, e.g. 4"
			read(*,*) tmparrint
			exit
		else if (ifuncsel==501.or.ifuncsel==502.or.ifuncsel==503.or.ifuncsel==504.or.ifuncsel==505) then !Information entropy related functions, need isolated atomic density
			call setpromol
			exit
		else
			if (ifuncsel==8) then
				inucespplot=1 !For deal with plotting nucesp property, special treatment is needed
			else
				inucespplot=0
			end if
			!Only for correlation hole/factor, source function, linear response kernel... reference point will be marked on contour map
			if (ifuncsel==17.or.ifuncsel==19.or.(ifuncsel==100.and.iuserfunc==24)) then
				imarkrefpos=1
			else
				imarkrefpos=0
			end if
			exit
		end if
	end do
end if

if (iorbsel2==0) then
    write(*,*) "-10 Return to main menu"
	write(*,*) "Draw which kind of map?"
	write(*,*) "1 Color-filled map (with/without contour lines)"
	write(*,*) "2 Contour line map"
	write(*,*) "3 Relief map"
	write(*,*) "4 Shaded relief map"
	write(*,*) "5 Shaded relief map with projection"
    if (itask==0) then
		write(*,*) "6 Gradient lines map with/without contour lines"
		write(*,*) "7 Vector field map with/without contour lines"
    end if
	read(*,*) idrawtype
else !Plotting two orbitals, only contour line map is available
	idrawtype=2
end if

if (idrawtype==1.or.idrawtype==2.or.idrawtype==6.or.idrawtype==7) then  !Initialize contour line setting
	if (idrawtype==2.or.idrawtype==6) idrawcontour=1
	if (idrawtype==1.or.idrawtype==7) idrawcontour=0
	if (idrawtype/=1) iatom_on_plane=1
	ilabel_on_contour=0
	if (ifuncsel==9.or.ifuncsel==10) then !A special contour setting suitable for ELF and LOL
        call gencontour(1,0D0,0D0,0)
	else !General contour setting for other real space functions
        call gencontour(0,0D0,0D0,0)
	end if
else if (idrawtype==-10) then
    return
end if

write(*,*) " -10 Return to main menu"
write(*,*) "How many grids in the two dimensions, respectively?"
if (itask==0.or.itask==2) then
	if (idrawtype==1.or.idrawtype==2.or.idrawtype==6) then
		if (ifdoESP(ifuncsel)) then
			write(*,*) "(100,100 is recommended)" !Because calculating ESP is very time consuming, so use lower grid
		else
			write(*,*) "(200,200 is recommended)"
		end if
	else if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
		write(*,*) "(100,100 is recommended)"
	else if (idrawtype==7) then
		write(*,*) "(80,80 is recommended)"
	end if
	write(*,*) "Hint: You can press ENTER button directly to use recommended value"
	read(*,"(a)") c200tmp

	if (c200tmp==' ') then !Pressing ENTER button directly
		if (idrawtype==1.or.idrawtype==2.or.idrawtype==6) then
			if (ifdoESP(ifuncsel)) then
				ngridnum1=100
				ngridnum2=100
			else
				ngridnum1=200
				ngridnum2=200
			end if
		else if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
			ngridnum1=100
			ngridnum2=100
		else if (idrawtype==7) then
			ngridnum1=80
			ngridnum2=80
		end if
	else if (c200tmp=='-10') then
		return
	else
		read(c200tmp,*) ngridnum1,ngridnum2
	end if
else
	write(*,*) "If press ENTER button directly, the recommended 100,100 will be used"
    read(*,"(a)") c200tmp
    if (c200tmp==" ") then
        ngridnum1=100
        ngridnum2=100
	else
		read(c200tmp,*) ngridnum1,ngridnum2
	end if
end if

if (allocated(planemat)) deallocate(planemat)
if (allocated(planemattmp)) deallocate(planemattmp)
allocate(planemat(ngridnum1,ngridnum2),planemattmp(ngridnum1,ngridnum2)) !planemattmp is used in many cases below
if (ncustommap/=0) then
	if (allocated(planemat_cust)) deallocate(planemat_cust)
	if (allocated(d1addtmp)) deallocate(d1addtmp)
	if (allocated(d1mintmp)) deallocate(d1mintmp)
	if (allocated(d2addtmp)) deallocate(d2addtmp)
	if (allocated(d2mintmp)) deallocate(d2mintmp)
	allocate(planemat_cust(ngridnum1,ngridnum2))
	allocate(d1addtmp(ngridnum1,ngridnum2))
	allocate(d1mintmp(ngridnum1,ngridnum2))
	allocate(d2addtmp(ngridnum1,ngridnum2))
	allocate(d2mintmp(ngridnum1,ngridnum2))
end if
if (idrawtype==6.or.idrawtype==7) then !Draw gradient lines
	if (allocated(d1add)) deallocate(d1add)
	if (allocated(d1min)) deallocate(d1min)
	if (allocated(d2add)) deallocate(d2add)
	if (allocated(d2min)) deallocate(d2min)
	if (allocated(gradd1)) deallocate(gradd1)
	if (allocated(gradd2)) deallocate(gradd2)
	allocate(d1add(ngridnum1,ngridnum2))
	allocate(d1min(ngridnum1,ngridnum2))
	allocate(d2add(ngridnum1,ngridnum2))
	allocate(d2min(ngridnum1,ngridnum2))
	allocate(gradd1(ngridnum1,ngridnum2))
	allocate(gradd2(ngridnum1,ngridnum2))
end if

orgx2D=0D0
orgy2D=0D0
orgz2D=0D0
v1x=0D0
v1y=0D0
v1z=0D0
v2x=0D0
v2y=0D0
v2z=0D0
transd1=0
transd2=0
rotplane=0
do while(.true.)
	write(*,*)
	write(*,*) " -10 Return to main menu" 
	write(*,*) "Please define the plane to be plotted"
	write(*,*) "1: XY plane   2: XZ plane   3: YZ plane"
	write(*,*) "4: Define by three atoms    5: Define by three given points"
	write(*,*) "6: Input origin and translation vector (For expert users)"
	write(*,*) "7: Parallel to a bond and meantime normal to a plane defined by three atoms"
	write(*,*) "8: Above or below the plane consisting of specific atoms"
	write(*,"(a,f8.4,a)") " 0: Set extension distance for plane type 1~5, current:",aug2D," Bohr"
	write(*,*) "-1: Set translation and rotation of the map for plane types 4, 5 and 8"
	read(*,*) plesel !Global variable
	aug2D2=aug2D !If don't draw gradient line map, needn't make the augment in the two dimension different
	if (plesel==-10) then
		return
	else if (plesel==0) then
		write(*,*) "Input extension distance in Bohr, e.g. 4.5"
        write(*,*) "Note: Negative value is also acceptable"
		read(*,*) aug2D
	else if (plesel==-1) then
		write(*,*) "Input translation of the content in X and Y directions, respectively (in Bohr)"
		write(*,*) "e.g. 0.5,-1.2"
        write(*,*) "If you press ENTER button directly, translation will not be applied"
        read(*,"(a)") c200tmp
        if (c200tmp==" ") then
			transd1=0
            transd2=0
        else
			read(c200tmp,*) transd1,transd2
        end if
		write(*,"(a)") " Input rotation angle of the plotting plane (in degree), positive and negative value correspond to clockwise and anticlockwise, respectively. e.g. 30.5"
		read(*,*) rotplane
		write(*,*) "Done!"
	else if (plesel==1) then
		write(*,*) "Input Z value in Bohr, e.g. 0.2"
		write(*,*) "Note: If the unit is in Angstrom, add ""a"" suffix, e.g. -1.6a"
		read(*,*) c200tmp
		if (index(c200tmp,'a')/=0) then
			read(c200tmp(1:len_trim(c200tmp)-1),*) orgz2D
			orgz2D=orgz2D/b2a
		else
			read(c200tmp,*) orgz2D
		end if
		if (ncenter==0) then
			orgx2D=orgx
			orgy2D=orgy
			v1x=(endx-orgx)/(ngridnum1-1)
			v2y=(endy-orgy)/(ngridnum2-1)
		else
			!Adjust aug2D/aug2D2 to make the length in two direction same
			!Because of bug in stream() of dislin 9.5D, if the two length unequal, some region won't be shown
			if (idrawtype==6.or.idrawtype==7) then  
				sup=(maxval(a%x)-minval(a%x))-(maxval(a%y)-minval(a%y))
				if (sup>0) aug2D2=aug2D+sup/2
				if (sup<0) aug2D=aug2D-sup/2
			end if
			orgx2D=minval(a%x)-aug2D
			orgy2D=minval(a%y)-aug2D2
			v1x=(maxval(a%x)+aug2D-orgx2D)/ngridnum1
			v2y=(maxval(a%y)+aug2D2-orgy2D)/ngridnum2
		end if
	else if (plesel==2) then
		write(*,*) "Input Y value in Bohr, e.g. 0.2"
		write(*,*) "Note: If the unit is in Angstrom, add ""a"" suffix, e.g. -1.6a"
		read(*,*) c200tmp
		if (index(c200tmp,'a')/=0) then
			read(c200tmp(1:len_trim(c200tmp)-1),*) orgy2D
			orgy2D=orgy2D/b2a
		else
			read(c200tmp,*) orgy2D
		end if
		if (ncenter==0) then
			orgx2D=orgx
			orgz2D=orgz
			v1x=(endx-orgx)/(ngridnum1-1)
			v2z=(endz-orgz)/(ngridnum2-1)
		else
			if (idrawtype==6.or.idrawtype==7) then !adjust aug2D/aug2D2
				sup=(maxval(a%x)-minval(a%x))-(maxval(a%z)-minval(a%z))
				if (sup>0) aug2D2=aug2D+sup/2
				if (sup<0) aug2D=aug2D-sup/2
			end if
			orgx2D=minval(a%x)-aug2D
			orgz2D=minval(a%z)-aug2D2
			v1x=(maxval(a%x)+aug2D-orgx2D)/ngridnum1
			v2z=(maxval(a%z)+aug2D2-orgz2D)/ngridnum2
		end if
	else if (plesel==3) then
		write(*,*) "Input X value in Bohr, e.g. 0.2"
		write(*,*) "Note: If the unit is in Angstrom, add ""a"" suffix, e.g. -1.6a"
		read(*,*) c200tmp
		if (index(c200tmp,'a')/=0) then
			read(c200tmp(1:len_trim(c200tmp)-1),*) orgx2D
			orgx2D=orgx2D/b2a
		else
			read(c200tmp,*) orgx2D
		end if
		if (ncenter==0) then
			orgy2D=orgy
			orgz2D=orgz
			v1y=(endy-orgy)/(ngridnum1-1)
			v2z=(endz-orgz)/(ngridnum2-1)
		else
			if (idrawtype==6.or.idrawtype==7) then !adjust aug2D/aug2D2
				sup=(maxval(a%y)-minval(a%y))-(maxval(a%z)-minval(a%z))
				if (sup>0) aug2D2=aug2D+sup/2
				if (sup<0) aug2D=aug2D-sup/2
			end if
			orgy2D=minval(a%y)-aug2D
			orgz2D=minval(a%z)-aug2D2
			v1y=(maxval(a%y)+aug2D-orgy2D)/ngridnum1
			v2z=(maxval(a%z)+aug2D2-orgz2D)/ngridnum2
		end if
	else if (plesel==4.or.plesel==5.or.plesel==8) then
		if (plesel==4) then
			write(*,*) "Input index of three atoms, e.g. 3,6,7"
			read(*,*) i1,i2,i3
			if (i1==i2.or.i1==i3.or.i2==i3.or.min(i1,i2,i3)<1.or.max(i1,i2,i3)>ncenter) then
				if (min(i1,i2,i3)<1.or.max(i1,i2,i3)>ncenter) write(*,*) "Atom indices are out of valid range, please input again"
				if (i1==i2.or.i1==i3.or.i2==i3) write(*,*) "Atom indices are duplicated, please input again"
				cycle
			end if
			a1x=a(i1)%x
			a1y=a(i1)%y
			a1z=a(i1)%z
			a2x=a(i2)%x
			a2y=a(i2)%y
			a2z=a(i2)%z
			a3x=a(i3)%x
			a3y=a(i3)%y
			a3z=a(i3)%z
		else if (plesel==5) then
			write(*,*) "Input x,y,z of point 1 (in Bohr, e.g. 1.2,1.3,0.0)"
            write(*,*) "To input coorindates in Angstrom, you can add "" A"" suffix, e.g. 0.23,0.5,-1 A"
            read(*,"(a)") c80tmp
			read(c80tmp,*) a1x,a1y,a1z
            if (index(c80tmp,'A')/=0) then
				a1x=a1x/b2a
				a1y=a1y/b2a
				a1z=a1z/b2a
            end if
			write(*,*) "Input x,y,z of point 2 (in Bohr)"
            write(*,*) "To input coorindates in Angstrom, you can add "" A"" suffix, e.g. 0.23,0.5,-1 A"
            read(*,"(a)") c80tmp
			read(c80tmp,*) a2x,a2y,a2z
            if (index(c80tmp,'A')/=0) then
				a2x=a2x/b2a
				a2y=a2y/b2a
				a2z=a2z/b2a
            end if
			write(*,*) "Input x,y,z of point 3 (in Bohr)"
            write(*,*) "To input coorindates in Angstrom, you can add "" A"" suffix, e.g. 0.23,0.5,-1 A"
            read(*,"(a)") c80tmp
			read(c80tmp,*) a3x,a3y,a3z
            if (index(c80tmp,'A')/=0) then
				a3x=a3x/b2a
				a3y=a3y/b2a
				a3z=a3z/b2a
            end if
		else if (plesel==8) then !This way is complicated, so use a separate subroutine, make it equivalent to defining three points
			call define_parallel_plane
            plesel=5
            aug2D=0
            aug2D2=0
		end if
		v1x=a1x-a2x
		v1y=a1y-a2y
		v1z=a1z-a2z
		v2x=a3x-a2x
		v2y=a3y-a2y
		v2z=a3z-a2z
		rnorm1=dsqrt(v1x**2+v1y**2+v1z**2) !Norm of vector 1
		rnorm2=dsqrt(v2x**2+v2y**2+v2z**2)
		if (abs(v1x*v2x+v1y*v2y+v1z*v2z)/(rnorm1*rnorm2)>0.999D0) then
			if (plesel==4) write(*,*) "The three atoms should not lie in the same line!"
			if (plesel==5) write(*,*) "The three points should not lie in the same line!"
			cycle
		end if
		rangle=acos( abs(v1x*v2x+v1y*v2y+v1z*v2z)/(rnorm1*rnorm2) )
        !For gradient line map plotted by "stream" of DISLIN, we must make lengths (d1,d2) of the two dimensions of the plotting
        !region identical here, so modify aug2D/aug2D2 to realize this aim. Because even if the passed arrays are completely correct and cover the whole
        !plotting region, gradient lines will still only appear in a cubic region whose side length is identical to min(dist1,dist2),
        !making a ugly blank region in the map. This should be a bug of DISLIN. I find it is impossible to overcome this by manually define
        !the so called starting point arrays in the "stream"
		if (idrawtype==6) then !adjust aug2D/aug2D2
			sup=rnorm1-rnorm2*cos(pi/2-rangle)
			if (sup>0) aug2D2=aug2D+sup/2
			if (sup<0) aug2D=aug2D-sup/2
		end if
		dist1=rnorm1+2*aug2D !Total length in direction 1
		dist2=rnorm2*cos(pi/2-rangle)+2*aug2D2 !Total length in direction 2
	!		write(*,*) "angle",acos(cos(pi/2-rangle))/pi*180
		d1=dist1/ngridnum1  !Transitional step length in direction 1
		d2=dist2/ngridnum2
		v1x=v1x*d1/rnorm1  !Make the norm of v1 equal to expected step lengh (d1)
		v1y=v1y*d1/rnorm1
		v1z=v1z*d1/rnorm1
		schmit=(v1x*v2x+v1y*v2y+v1z*v2z)/(v1x**2+v1y**2+v1z**2) !Use schmit method to make v2 ortho to v1
		v2x=v2x-schmit*v1x
		v2y=v2y-schmit*v1y
		v2z=v2z-schmit*v1z
		rnorm2=dsqrt(v2x**2+v2y**2+v2z**2)
		v2x=v2x*d2/rnorm2   !Make the norm of v2 equal to expected step lengh (d2)
		v2y=v2y*d2/rnorm2
		v2z=v2z*d2/rnorm2
		orgx2D=a2x-aug2D/d1*v1x-aug2D2/d2*v2x  !aug2D/d1*v1x=aug2D*(v1x/d1), v1x/d1 correspond the x component of unit vector in v1x direction
		orgy2D=a2y-aug2D/d1*v1y-aug2D2/d2*v2y
		orgz2D=a2z-aug2D/d1*v1z-aug2D2/d2*v2z
		
		!Customized translation. Currently d1 and d2 are norm of the two translation vectors
		if (transd1/=0.or.transd2/=0) then
			orgx2D=orgx2D-transd1*v1x/d1-transd2*v2x/d2
			orgy2D=orgy2D-transd1*v1y/d1-transd2*v2y/d2
			orgz2D=orgz2D-transd1*v1z/d1-transd2*v2z/d2
		end if
		!Customized rotation
		if (rotplane/=0) then
			xcen=orgx2D+dist1/2D0*v1x/d1+dist2/2D0*v2x/d2
			ycen=orgy2D+dist1/2D0*v1y/d1+dist2/2D0*v2y/d2
			zcen=orgz2D+dist1/2D0*v1z/d1+dist2/2D0*v2z/d2
			vec1old(1)=v1x;vec1old(2)=v1y;vec1old(3)=v1z
			vec2old(1)=v2x;vec2old(2)=v2y;vec2old(3)=v2z
			vec1old=vec1old/dsqrt(sum(vec1old**2)) !Normalize
			vec2old=vec2old/dsqrt(sum(vec2old**2))
			rotplane=-rotplane/180*pi !If the graph is rotated clockwisely, the axis should be rotated anticlockwisely
			vec1=cos(rotplane)*vec1old-sin(rotplane)*vec2old !Matrix rotation
			vec2=cos(rotplane)*vec2old+sin(rotplane)*vec1old
			vec1=vec1*d1
			vec2=vec2*d2
			v1x=vec1(1);v1y=vec1(2);v1z=vec1(3)
			v2x=vec2(1);v2y=vec2(2);v2z=vec2(3)
			!Make center of the original plane map corresponds to the center of the new plane map
			orgx2D=xcen-dist1/2D0*v1x/d1-dist2/2D0*v2x/d2
			orgy2D=ycen-dist1/2D0*v1y/d1-dist2/2D0*v2y/d2
			orgz2D=zcen-dist1/2D0*v1z/d1-dist2/2D0*v2z/d2
		end if
	else if (plesel==6) then
		write(*,*) "Input origin of x,y,z in Bohr, e.g. 3.5,-1,0.2"
		read(*,*) orgx2D,orgy2D,orgz2D
		write(*,*) "Input x,y,z of transitional vector 1 in Bohr, e.g. 0.08,0.03,0"
		read(*,*) v1x,v1y,v1z
		write(*,*) "Input x,y,z of transitional vector 2 in Bohr, should be orthogonal to vector 1"
		read(*,*) v2x,v2y,v2z
		d1=dsqrt(v1x**2+v1y**2+v1z**2)
		d2=dsqrt(v2x**2+v2y**2+v2z**2)
		dist1=d1*ngridnum1
		dist2=d2*ngridnum2
		a1x=orgx2D !Although a1x...a3z is no use for generate data, but these are critical for plotting atom label (subroutine drawplane)
		a1y=orgy2D
		a1z=orgz2D
		a2x=orgx2D+ngridnum1*v1x
		a2y=orgy2D+ngridnum1*v1y
		a2z=orgz2D+ngridnum1*v1z
		a3x=orgx2D+ngridnum2*v2x
		a3y=orgy2D+ngridnum2*v2y
		a3z=orgz2D+ngridnum2*v2z
	else if (plesel==7) then
		write(*,*) "Input two atoms to define the bond, e.g. 4,5"
		read(*,*) iatm1,iatm2
		write(*,*) "Input three atoms to define a plane, e.g. 1,4,7"
		read(*,*) jatm1,jatm2,jatm3
		write(*,*) "Input length of X-axis in Bohr, e.g. 10"
		read(*,*) dist1
		write(*,*) "Input length of Y-axis in Bohr, e.g. 8"
		read(*,*) dist2
		v1x=a(iatm2)%x-a(iatm1)%x
		v1y=a(iatm2)%y-a(iatm1)%y
		v1z=a(iatm2)%z-a(iatm1)%z
		call pointABCD(a(jatm1)%x,a(jatm1)%y,a(jatm1)%z,a(jatm2)%x,a(jatm2)%y,a(jatm2)%z,a(jatm3)%x,a(jatm3)%y,a(jatm3)%z,v2x,v2y,v2z,tmpD)
		schmit=(v1x*v2x+v1y*v2y+v1z*v2z)/(v1x**2+v1y**2+v1z**2) !Use schmit method to make v2 ortho to v1
		v2x=v2x-schmit*v1x
		v2y=v2y-schmit*v1y
		v2z=v2z-schmit*v1z
		d1=dist1/(ngridnum1-1)
		d2=dist2/(ngridnum2-1)
		rnorm1=dsqrt(v1x**2+v1y**2+v1z**2)
		v1x=v1x/rnorm1*d1
		v1y=v1y/rnorm1*d1
		v1z=v1z/rnorm1*d1
		rnorm2=dsqrt(v2x**2+v2y**2+v2z**2)
		v2x=v2x/rnorm2*d2
		v2y=v2y/rnorm2*d2
		v2z=v2z/rnorm2*d2
		orgx2D=(a(iatm1)%x+a(iatm2)%x)/2 - v1x*ngridnum1/2 - v2x*ngridnum2/2
		orgy2D=(a(iatm1)%y+a(iatm2)%y)/2 - v1y*ngridnum1/2 - v2y*ngridnum2/2
		orgz2D=(a(iatm1)%z+a(iatm2)%z)/2 - v1z*ngridnum1/2 - v2z*ngridnum2/2
		a1x=orgx2D !Although a1x...a3z is no use for generate data, but these are critical for plotting atom label (subroutine drawplane)
		a1y=orgy2D
		a1z=orgz2D
		a2x=orgx2D+ngridnum1*v1x
		a2y=orgy2D+ngridnum1*v1y
		a2z=orgz2D+ngridnum1*v1z
		a3x=orgx2D+ngridnum2*v2x
		a3y=orgy2D+ngridnum2*v2y
		a3z=orgz2D+ngridnum2*v2z
	end if
	
	if (plesel/=0.and.plesel/=-1) exit
end do

write(*,*)
write(*,"(' X/Y/Z of origin of the plane: ',3f10.5,' Bohr')") orgx2D,orgy2D,orgz2D
endx2D=orgx2D+v1x*(ngridnum1-1)+v2x*(ngridnum2-1)
endy2D=orgy2D+v1y*(ngridnum1-1)+v2y*(ngridnum2-1)
endz2D=orgz2D+v1z*(ngridnum1-1)+v2z*(ngridnum2-1)
rnorm1=dsqrt(v1x**2+v1y**2+v1z**2) !The final length of vector 1
rnorm2=dsqrt(v2x**2+v2y**2+v2z**2) !The final length of vector 2
write(*,"(' X/Y/Z of end of the plane:    ',3f10.5,' Bohr')") endx2D,endy2D,endz2D
write(*,"(' X/Y/Z of translation vector 1:',3f9.5,' Bohr, Norm:',f9.5)") v1x,v1y,v1z,rnorm1
write(*,"(' X/Y/Z of translation vector 2:',3f9.5,' Bohr, Norm:',f9.5)") v2x,v2y,v2z,rnorm2

!The infinitesimal in each direction for gradient plot
diff=1D-5
diffv1x=diff*v1x/rnorm1
diffv1y=diff*v1y/rnorm1
diffv1z=diff*v1z/rnorm1
diffv2x=diff*v2x/rnorm2
diffv2y=diff*v2y/rnorm2
diffv2z=diff*v2z/rnorm2

if (itask==1) then !NICS-2D scan
	do while(.true.)
		write(*,*)
		write(*,*) "1 Generate Gaussian input file for NICS-2D scanning"
		write(*,*) "2 Load Gaussian output file of NICS-2D scanning"
		read(*,*) isel
		if (isel==1) then
			write(*,*) "Input the path of Gaussian template file of performing NMR task"
			write(*,*) "e.g. D:\Aqours\Mari\shiny.gjf"
			write(*,"(a)") " Note: In this file, the coordinate part should be recorded as [geometry], which will be automatically replaced with current geometry"
			do while(.true.)
				read(*,"(a)") c2000tmp
				inquire(file=c2000tmp,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again!"
			end do

			open(10,file=c2000tmp,status="old")
			open(11,file="NICS_2D.gjf",status="replace")
			do while (.true.)
				read(10,"(a)",iostat=ierror) c2000tmp
				if (ierror/=0) exit
				if (index(c2000tmp,"#")/=0) then
					if (index(c2000tmp,"geom=conn")==0) then
						write(11,"(a)") trim(c2000tmp)//" geom=connectivity"
					else
						write(11,"(a)") trim(c2000tmp)
					end if
				else if (index(c2000tmp,"[geometry]")/=0) then
					do iatm=1,ncenter
						write(11,"(a,3f12.6)") ind2name(a(iatm)%index),a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
					end do
					do ipt=1,ngridnum1
						do jpt=1,ngridnum2
							call get2Dgridxyz(ipt,jpt,rnowx,rnowy,rnowz)
							write(11,"('Bq',3f12.6)") rnowx*b2a,rnowy*b2a,rnowz*b2a
						end do
					end do
					write(11,*)
					do itmp=1,ncenter+ngridnum1*ngridnum2
						write(11,"(i6)") itmp
					end do
				else
					write(11,"(a)") trim(c2000tmp)
				end if
			end do
			write(11,*)
			write(11,*)
			close(10)
			close(11)
			write(*,"(a)") " NICS_2D.gjf has been generated in current folder! Please check it, and then run it by Gaussian manually"
		else if (isel==2) then
			write(*,*) "Input path of Gaussian output file of NICS-2D scanning task"
			write(*,*) "e.g. D:\Aqours\Mari\shiny.out"
			do while(.true.)
				read(*,"(a)") c2000tmp
				inquire(file=c2000tmp,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again!"
			end do
			exit
		end if
    end do
	!Load shielding tensor of ghost atoms
	write(*,*) "Obtain which kind of NICS?"
    write(*,*) "0: Projection along specific direction"
	write(*,*) "1: Isotropic  2: Anisotropy  3: XX component  4: YY component  5: ZZ component"
    read(*,*) iload
    if (iload==0) then
		write(*,*) "Input direction vector, e.g. 3.2,1.1,-1.9"
        read(*,*) tmpvec
        tmpvec=tmpvec/dsqrt(sum(tmpvec**2))
    end if
	open(10,file=c2000tmp,status="old")
	call loclabel(10,"Isotropic =",ifound)
	if (ifound==0) then
		close(10)
		write(*,"(a)") " Error: Unable to find magnetic shielding tensor in this file! Please check keywords. Press ENTER button to return"
		read(*,*)
		return
	end if
	call loclabel(10,"Bq",ifound,0)
	do ipt=1,ngridnum1
		do jpt=1,ngridnum2
			read(10,"(a)") c80tmp
			read(c80tmp(26:),*) tiso
			read(c80tmp(52:),*) taniso
			read(10,*) c80tmp,tmpmat(1,1),c80tmp,tmpmat(1,2),c80tmp,tmpmat(1,3)
			read(10,*) c80tmp,tmpmat(2,1),c80tmp,tmpmat(2,2),c80tmp,tmpmat(2,3)
			read(10,*) c80tmp,tmpmat(3,1),c80tmp,tmpmat(3,2),c80tmp,tmpmat(3,3)
			read(10,*)
            if (iload==0) then
				planemat(ipt,jpt)=-prjmat(tmpmat(:,:),tmpvec(:))
            else if (iload==1) then
				planemat(ipt,jpt)=-tiso
			else if (iload==2) then
				planemat(ipt,jpt)=-taniso
			else if (iload==3) then
				planemat(ipt,jpt)=-tmpmat(1,1)
			else if (iload==4) then
				planemat(ipt,jpt)=-tmpmat(2,2)
			else if (iload==5) then
				planemat(ipt,jpt)=-tmpmat(3,3)
            end if
        end do
	end do
	close(10)
	write(*,*) "Loading finished!"

else if (itask==2) then !Hyper(polarizability) density
	if (abs(itype)==1) ntime=2
	if (abs(itype)==2) ntime=3
	if (abs(itype)==3) ntime=4
	do itime=1,ntime
        call dealloall(0)
        write(*,*) "Calculating "//trim(custommapname(itime))
        call readinfile(custommapname(itime),1)
		do ipt=1,ngridnum1
			do jpt=1,ngridnum2
				call get2Dgridxyz(ipt,jpt,tmpx,tmpy,tmpz)
                planemattmp(ipt,jpt)=fdens(tmpx,tmpy,tmpz)
			end do
		end do
        if (abs(itype)==1) then !Polarizability density
            if (itime==1) then
                planemat=-planemattmp
            else if (itime==2) then
                planemat=planemat+planemattmp
            end if
        else if (abs(itype)==2) then !First hyperpolarizability density
            if (itime==1) then
                planemat=planemattmp
            else if (itime==2) then
                planemat=planemat-2*planemattmp
            else if (itime==3) then
                planemat=planemat+planemattmp
            end if
        else if (abs(itype)==3) then !Second hyperpolarizability density
            if (itime==1) then
                planemat=-planemattmp
            else if (itime==2) then
                planemat=planemat+2*planemattmp
            else if (itime==3) then
                planemat=planemat-2*planemattmp
            else if (itime==4) then
                planemat=planemat+planemattmp
            end if
        end if
    end do
    fieldstr=0.003D0
    if (abs(itype)==1) then
        planemat=planemat/(2*fieldstr)
    else if (abs(itype)==2) then
        planemat=planemat/fieldstr**2
    else if (abs(itype)==3) then
        planemat=planemat/(2*fieldstr**3)
    end if
    !Transform (hyper)polarizability to spatial contribution
    if (itype<0) then
		do ipt=1,ngridnum1
			do jpt=1,ngridnum2
				call get2Dgridxyz(ipt,jpt,tmpx,tmpy,tmpz)
                if (itype2==1) then
                    planemat(ipt,jpt)=-tmpx*planemat(ipt,jpt)
                else if (itype2==2) then
                    planemat(ipt,jpt)=-tmpy*planemat(ipt,jpt)
                else if (itype2==3) then
                    planemat(ipt,jpt)=-tmpz*planemat(ipt,jpt)
                end if
			end do
		end do
    end if
    call dealloall(0)
    write(*,*) "Reloading "//trim(firstfilename)
    call readinfile(firstfilename,1)
	
else if (iplaneextdata==1) then !Export plane data to external file, and then load data from it
	open(10,file="planept.txt",status="replace")
	open(11,file="cubegenpt.txt",status="replace")
	write(10,*) ngridnum1*ngridnum2
	do ipt=1,ngridnum1
		do jpt=1,ngridnum2
			call get2Dgridxyz(ipt,jpt,rnowx,rnowy,rnowz)
			write(10,"(3f16.8)") rnowx,rnowy,rnowz
			write(11,"(3f16.8)") rnowx*b2a,rnowy*b2a,rnowz*b2a
		end do
	end do
	close(10)
	close(11)
	write(*,"(/,a)") " The coordinate of all points needed to be calculated have been outputted to plane.txt in current folder, the unit is in Bohr"
	write(*,"(a)") " cubegenpt.txt is also outputted, which is similar to plane.txt, but the unit is in Angstrom, &
	&and there is no first line (the number of points). It can be directly utilized by cubegen"
	write(*,"(a)") " For example ""cubegen 0 potential CNT.fch result.cub -5 h < cubegenpt.txt"""
	write(*,*)
	write(*,"(a)") " Now input the path of the file containing function values, e.g. C:\t.txt, whose format should be identical to plane.txt, but with function value in the fourth column"
	write(*,"(a)") " Note: If the suffix is .cub, then the file will be recognized and loaded as output file of cubegen"
	do while(.true.)
		read(*,"(a)") c200tmp
		inquire(file=c200tmp,exist=alive)
		if (alive) exit
		write(*,*) "File cannot be found, input again"
	end do
	open(10,file=c200tmp,status="old")
	if (index(c200tmp,".cub")/=0) then !cubegen output
		do iskip=1,6+ncenter
			read(10,*)
		end do
	else
		read(10,*)
	end if
	do ipt=1,ngridnum1
		do jpt=1,ngridnum2
			read(10,*) rnouse,rnouse,rnouse,planemat(ipt,jpt)
		end do
	end do
	close(10)
	
else !Common case, start calculation of plane data
	write(*,*)
	if (ifuncsel/=4) call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation, but don't do this when analyzing MO
	call gen_GTFuniq(1) !Generate unique GTFs, for faster evaluation in orbderv
	write(*,*) "Calculating plane data, please wait..."	
	if (ifuncsel/=12) then
        if (ifPBC>0) then
	        if (expcutoff_PBC<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' will be ignored')") expcutoff_PBC
        else
    	    if (expcutoff<0) write(*,"(' Note: All exponential functions exp(x) with x<',f8.3,' will be ignored')") expcutoff
        end if
    end if
	call walltime(iwalltime1)
	icustom=0
	planemat=0D0 !For promolecular property, first clean up planemat
	if (ipromol==1) goto 401 ! To obtain promolecule property, pass the first loaded molecule
	!Note: If the task refers to plotting ESP map, the function of drawing gradient line will be disabled

	400	if (ifuncsel==12.and.idrawtype/=6.and.idrawtype/=7) then
		call planeesp !Special treatment to calculate ESP
	else if (ifuncsel==112) then !Calculate atomic Hirshfeld weighting function
		call genhirshplanewei(tmparrint,size(tmparrint),iHirshdenstype)
		ncustommap=0
	else if (ifuncsel==500.or.ifuncsel==510.or.ifuncsel==511.or.ifuncsel==512) then
	!500: Calculate rho(A)*ln[rho(A)/rho0(A)], 510: Calculate rho(A), 511: Calculate rho0(A), 512: other
		call genhirshplanewei(tmparrint,size(tmparrint),1)
		ncustommap=0
		do i=1,ngridnum1 !Now planemat is Hirshfeld weight of iatmentropy, and planemattmp is its density in free-state
			do j=1,ngridnum2
				call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
				rhoA=planemat(i,j)*fdens(rnowx,rnowy,rnowz)
				if (ifuncsel==500) planemat(i,j)=rhoA*log(rhoA/planemattmp(i,j))
				if (ifuncsel==510) planemat(i,j)=rhoA
				if (ifuncsel==511) planemat(i,j)=planemattmp(i,j)
				if (ifuncsel==512) planemat(i,j)=log(rhoA/planemattmp(i,j))
			end do
		end do
	else if (ifuncsel==501) then !Calculate sum{rho(A)*ln[rho(A)/rho0(A)]}
		call genentroplane(1)
		ncustommap=0
	else if (ifuncsel==502) then !Calculate sum(x), where x=[rho(A)-rho0(A)]/rho(A)
		call genentroplane(2)
		ncustommap=0
	else if (ifuncsel==503) then !Calculate difference between total relative Shannon entropy and deformation density 
		call genentroplane(3)
		ncustommap=0
	else if (ifuncsel==504) then !Calculate 2nd relative Onicescu information sum{[rho(A)]^2/rho0(A)}
		call genentroplane(4)
		ncustommap=0
	else if (ifuncsel==505) then !Calculate 3rd relative Onicescu information sum{[rho(A)]^3/[rho0(A)]^2}/2
		call genentroplane(5)
		ncustommap=0
    else if (ifuncsel==100.and.(iuserfunc==57.or.iuserfunc==58.or.iuserfunc==59)) then !Calculate g1,g2,g3 terms defined by Shubin, they rely on rho_0
        call g1g2g3plane
		ncustommap=0
	else !Common case
        nthreads_old=nthreads
        !For some ESP related functions, initialize LIBRETA so that fast code will be used
        if (ifdoESP(ifuncsel).and.(iESPcode==2.or.iESPcode==3)) then
            call doinitlibreta(1)
            if (isys==1.and.nthreads>12) nthreads=12
        end if
        iprog=0
	    !$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz) shared(iprog,planemat,d1add,d1min,d2add,d2min) schedule(dynamic) NUM_THREADS(nthreads)
		do i=1,ngridnum1
			do j=1,ngridnum2
				call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
				if (ifuncsel==111) then
					planemat(i,j)=beckewei(rnowx,rnowy,rnowz,iatmbecke1,iatmbecke2)
				else
					planemat(i,j)=calcfuncall(ifuncsel,rnowx,rnowy,rnowz)
					if (ifuncsel==4.and.iorbsel2/=0) planemattmp(i,j)=fmo(rnowx,rnowy,rnowz,iorbsel2) !Calculate another orbital together
				end if
				if (idrawtype==6.or.idrawtype==7) then !Generate two vector to plot vector field line graph
					d1add(i,j)=calcfuncall(ifuncsel,rnowx+diffv1x,rnowy+diffv1y,rnowz+diffv1z)
					d1min(i,j)=calcfuncall(ifuncsel,rnowx-diffv1x,rnowy-diffv1y,rnowz-diffv1z)
					d2add(i,j)=calcfuncall(ifuncsel,rnowx+diffv2x,rnowy+diffv2y,rnowz+diffv2z)
					d2min(i,j)=calcfuncall(ifuncsel,rnowx-diffv2x,rnowy-diffv2y,rnowz-diffv2z)
				end if
			end do
			!$OMP CRITICAL
            iprog=iprog+1
			call showprog(iprog,ngridnum1)
            !$OMP END CRITICAL
		end do
	    !$OMP END PARALLEL DO
        nthreads=nthreads_old
	end if

	401	if (ncustommap/=0) then !Calculate data for custom map
		if (icustom==0) then !The first time
			planemat_cust=planemat
			if (idrawtype==6.or.idrawtype==7) then
				d1addtmp=d1add
				d1mintmp=d1min
				d2addtmp=d2add
				d2mintmp=d2min
			end if
		else if (icustom/=0) then
			if (customop(icustom)=='+') then
				planemat_cust=planemat_cust+planemat
				if (idrawtype==6.or.idrawtype==7) then
					d1addtmp=d1addtmp+d1add
					d1mintmp=d1mintmp+d1min
					d2addtmp=d2addtmp+d2add
					d2mintmp=d2mintmp+d2min
				end if
			else if (customop(icustom)=='-') then
				planemat_cust=planemat_cust-planemat
				if (idrawtype==6.or.idrawtype==7) then
					d1addtmp=d1addtmp-d1add
					d1mintmp=d1mintmp-d1min
					d2addtmp=d2addtmp-d2add
					d2mintmp=d2mintmp-d2min
				end if
			else if (customop(icustom)=='x'.or.customop(icustom)=='*') then
				planemat_cust=planemat_cust*planemat
				if (idrawtype==6.or.idrawtype==7) then
					d1addtmp=d1addtmp*d1add
					d1mintmp=d1mintmp*d1min
					d2addtmp=d2addtmp*d2add
					d2mintmp=d2mintmp*d2min
				end if
			else if (customop(icustom)=='/') then
				planemat_cust=planemat_cust/planemat
				if (idrawtype==6.or.idrawtype==7) then
					d1addtmp=d1addtmp/d1add
					d1mintmp=d1mintmp/d1min
					d2addtmp=d2addtmp/d2add
					d2mintmp=d2mintmp/d2min
				end if
			end if
		end if
		if (icustom/=ncustommap) then !Not the final time
			icustom=icustom+1
			filename=custommapname(icustom)
            call savePBCinfo
			call dealloall(0)
			write(*,"(' Calculating: ',a)") trim(filename)
			call readinfile(filename,1)
            call loadPBCinfo
			if (ifuncsel/=4) call delvirorb(0)
			call gen_GTFuniq(1) !Generate unique GTFs, for faster evaluation in orbderv
			!Input the MO index for current file. Since the MO index may be not the same as the first loaded one
			if (ifuncsel==4) then
				write(*,"(' Input index of the orbital to be calculated for ',a,', e.g. 3')") trim(filename)
				read(*,*) iorbsel
			end if
			goto 400
		else !The final time, reload the first loaded system
			planemat=planemat_cust
			if (idrawtype==6.or.idrawtype==7) then
				d1add=d1addtmp
				d1min=d1mintmp
				d2add=d2addtmp
				d2min=d2mintmp
			end if
			call dealloall(0)
			write(*,"(' Reloading: ',a)") trim(firstfilename)
			call readinfile(firstfilename,1)
		end if
	end if

	if (idrawtype==6.or.idrawtype==7) then !Finish the finite differential to yield gradient
		gradd1=(d1add-d1min)/2/diff
		gradd2=(d2add-d2min)/2/diff
	end if
	
	call del_GTFuniq !Destory unique GTF informtaion
    call delvirorb_back(1) !delvirorb may have taken effect, now restore to previous wavefunction
	call walltime(iwalltime2)
	write(*,"(/,' Calculation took up wall clock time',i10,' s',/)") iwalltime2-iwalltime1
end if

valmin=minval(planemat)
valmax=maxval(planemat)
write(*,*) "The minimum of data:",valmin
write(*,*) "The maximum of data:",valmax

! Set default lower and upper of color scale (Z-axis) for plot
surcolorzmin=-3
surcolorzmax=3
drawlowlim=valmin
drawuplim=valmax
if (itask==1) then !NICS-2D
	drawlowlim=-30
    drawuplim=30
else if (itask==2) then !NICS-2D
	drawlowlim=-10
    drawuplim=10
else if (ifuncsel==1) then
	drawlowlim=0D0
	drawuplim=0.65D0
else if (ifuncsel==2) then
	drawlowlim=0D0
	drawuplim=0.65D0
else if (ifuncsel==3) then 
	drawlowlim=-8D0
	drawuplim=15D0
else if (ifuncsel==4) then
	drawlowlim=-0.8D0
	drawuplim=0.8D0
else if (ifuncsel==5) then
	drawlowlim=-0.1D0
	drawuplim=0.1D0
else if (ifuncsel==8) then !Nuclear ESP
	if (ifiletype==4) then !Atomic charge ESP
		drawlowlim=-0.4D0
		drawuplim=0.4D0
	else
		drawlowlim=0D0
		drawuplim=50D0
		surcolorzmin=0
		surcolorzmax=50
	end if
else if (ifuncsel==9) then
	drawlowlim=0D0
	drawuplim=1D0
else if (ifuncsel==10) then
	drawlowlim=0D0
	drawuplim=0.8D0
else if (ifuncsel==11) then
	drawlowlim=0D0
	drawuplim=0.1D0
else if (ifuncsel==12) then
	drawlowlim=-0.1D0
	drawuplim=0.1D0
else if (ifuncsel==13.or.ifuncsel==14) then
	drawlowlim=0D0
	drawuplim=1D0
else if (ifuncsel==15.or.ifuncsel==16) then
	drawlowlim=-0.65D0
	drawuplim=0.65D0
else if (ifuncsel==17) then
	drawlowlim=-0.5D0
	drawuplim=0.1D0
else if (ifuncsel==18) then
	drawlowlim=0D0
	drawuplim=2D0
else if (ifuncsel==22.or.ifuncsel==23) then !delta-g
	drawlowlim=0D0
	drawuplim=0.5D0
else if (ifuncsel==24) then !IRI
	drawlowlim=0D0
	drawuplim=2D0
else if (ifuncsel==25) then !vdW potential
	drawlowlim=-1D0
	drawuplim=1D0
else if (ifuncsel==44) then !Orbital probability density
	drawlowlim=-1D-3
	drawuplim=0.05D0
else if (ifuncsel==111.or.ifuncsel==112) then !Becke/Hirshfeld weight
	drawlowlim=0D0
	drawuplim=1D0
else if (ifuncsel==100) then
	if (iuserfunc==20.or.iuserfunc==99) then !DORI/IRI
		drawlowlim=0D0
		drawuplim=2D0
	else if (iuserfunc>=910.and.iuserfunc<=914) then !Atomic weighting functions
		drawlowlim=0D0
		drawuplim=1.000001D0
    end if
end if
!Set up range of X and Y axes
if (plesel==1) then
	axlow1=orgx2D
	axhigh1=endx2D
	axlow2=orgy2D
	axhigh2=endy2D
else if (plesel==2) then
	axlow1=orgx2D
	axhigh1=endx2D
	axlow2=orgz2D
	axhigh2=endz2D
else if (plesel==3) then
	axlow1=orgy2D
	axhigh1=endy2D
	axlow2=orgz2D
	axhigh2=endz2D
else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
	axlow1=0D0
	axhigh1=dist1-d1
	axlow2=0D0
	axhigh2=dist2-d2
end if
!Step size between labels
planestpx=(axhigh1-axlow1)/7
planestpy=(axhigh2-axlow2)/7
planestpz=(drawuplim-drawlowlim)/10

XVU=150D0 !Reinitialize view
YVU=30D0
ZVU=7D0 !More suitable than 6D0 for drawing plane
iatmcontri=0 !=0 means haven't define atomic contribution
idrawintbasple=0 !Refresh (3,-1) information
nple3n1path=0
cp2ple3n1path=0
if (allocated(ple3n1path)) deallocate(ple3n1path)
i=-1

if (numcp>0) then
	write(*,"(/,a)") " If you did not modify ""CP_RGB_2D"" in settings.ini, then brown, blue, orange and green dots in this map (if any) &
    &correspond to position of (3,-3), (3,-1), (3,+1) and (3,+3) type of critical points, respectively"
end if

!----------------- post-processing menu, use plane data to draw graph
do while(.true.)
!! Because the max cycle is ngridnum1/2 - 1, so the upper coordinate likes dist1-d1 rather than dist1
	if ((i==-1.and.isilent==0).or.isavepic==1) then
		if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then !Draw 3D plane, first use drawplaneGUI to setup GUI, which then invokes drawplane
			call drawplanegui(axlow1,axhigh1,axlow2,axhigh2,drawlowlim,drawuplim,idrawtype)
		else
			call drawplane(axlow1,axhigh1,axlow2,axhigh2,drawlowlim,drawuplim)
		end if
		if (isavepic==1) write(*,"(a,a,a)") " Graph have been saved as ",trim(graphformat)," format with ""dislin"" prefix in current directory"
		isavepic=0
	end if
	
	!! After show the plot once, ask user what to do next. 0 and negative options are general options
	write(*,*)
    call menutitle("Post-processing menu of plotting plane map",8,1)
! 	write(*,*) "-10 Multiply data by the data in a plane text file"
	if (iatmcontri==0) write(*,*) "-9 Only plot the data around certain atoms"
	if (iatmcontri==1) write(*,*) "-9 Recover original plane data"
	if (ilenunit2D==1) write(*,*) "-8 Change length unit of the graph to Angstrom"
	if (ilenunit2D==2) write(*,*) "-8 Change length unit of the graph to Bohr"		
	write(*,*) "-7 Multiply the current data by a factor"
	write(*,*) "-6 Export the current plane data to plane.txt in current folder"
	write(*,*) "-5 Return to main menu"
    write(*,*) "-4 Save (load) all plotting settings to (from) an external file"
	write(*,*) "-3 Change other plotting settings"
	if (idrawtype==1.or.idrawtype==4.or.idrawtype==5) then
		write(*,"(a,2f7.3,f12.6)") " -2 Set label intervals in X, Y and color scale axes, current:",planestpx,planestpy,planestpz
	else
		write(*,"(a,2f7.3)") " -2 Set label intervals in X and Y axes, current:",planestpx,planestpy
	end if
	write(*,*) "-1 Show the graph again"
	write(*,*) "0 Save the graph to a graphical file in current folder"

	if (idrawtype==1) then !Color-filled map
		if (abs(drawlowlim)<1000000.and.abs(drawuplim)<1000000) then
			write(*,"(a,2f15.7)") " 1 Set lower&upper limit of color scale, current:",drawlowlim,drawuplim
		else
			write(*,"(a,2(1PE15.6))") " 1 Set lower&upper limit of color scale, current:",drawlowlim,drawuplim
		end if
		if (idrawcontour==1) write(*,*) "2 Disable showing contour lines"
		if (idrawcontour==0) write(*,*) "2 Enable showing contour lines"
		write(*,*) "3 Change contour line setting"
		if (iatom_on_plane==0) write(*,*) "4 Enable showing atom labels and reference point"
		if (iatom_on_plane==1) write(*,*) "4 Disable showing atom labels and reference point"
		if (numcp>0.or.numpath>0) write(*,*) "5 Set details of plotting critical points and paths"
		if (ibond_on_plane==0) write(*,*) "8 Enable showing bonds"
		if (ibond_on_plane==1) write(*,*) "8 Disable showing bonds"
		if (idrawplanevdwctr==0.and.iorbsel2==0.and.allocated(b)) write(*,*) "15 Enable showing the contour line corresponding to vdW surface (rho=0.001)" !meaningless if custom operation is performed
		if (idrawplanevdwctr==1) write(*,*) "15 Disable showing the contour line corresponding to vdW surface (rho=0.001)"
		if (idrawplanevdwctr==1) write(*,*) "16 Set label size, style and color of the contour line of vdW surface" !When iorbsel2/=0, that means plot another orbital, cubmattmp will be pre-occupied
		if (iatom_on_plane==1) then
			if (idrawintbasple==0) then
				write(*,"(a,f7.3,' Bohr')") " 17 Set distance threshold for showing atom labels, current:",disshowlabel
			else
				write(*,"(a,f6.3,' Bohr')") " 17 Set distance threshold for showing atoms & interbasin paths:",disshowlabel
			end if
		end if
		if (iatmlabtype==1) write(*,*) "18 Change style of atomic labels: Only plot element symbol"
		if (iatmlabtype==2) write(*,*) "18 Change style of atomic labels: Only plot atomic index"
		if (iatmlabtype==3) write(*,*) "18 Change style of atomic labels: Plot both element symbol and atomic index"
        write(*,"(a,a)") " 19 Set color transition, current: ",trim(clrtransname(iclrtrans))
	else if (idrawtype==2.or.idrawtype==6.or.idrawtype==7) then !Contour map, gradient line, vector field with/without contour
		if (iatom_on_plane==0) write(*,*) "1 Enable showing atom labels and reference point"
		if (iatom_on_plane==1) write(*,*) "1 Disable showing atom labels and reference point"
		if (ilabel_on_contour==0) write(*,*) "2 Enable showing isovalue on contour lines"
		if (ilabel_on_contour==1) write(*,*) "2 Disable showing isovalue on contour lines"
		write(*,*) "3 Change setting of contour lines"
		if (numcp>0.or.numpath>0) write(*,*) "4 Set details of plotting critical points and paths"
		if (idrawcontour==1) write(*,*) "5 Disable showing contour lines"
		if (idrawcontour==0) write(*,*) "5 Enable showing contour lines"
		if (numcp>0.and.idrawintbasple==0) write(*,*) "6 Generate and show interbasin paths"
		if (numcp>0.and.idrawintbasple==1) write(*,*) "6 Delete interbasin paths"
		if (numcp>0.and.idrawintbasple==0) write(*,*) "7 Set stepsize and maximal iteration for interbasin path generation"
		if (ibond_on_plane==0) write(*,*) "8 Enable showing bonds"
		if (ibond_on_plane==1) write(*,*) "8 Disable showing bonds"
        if (idrawcontour==1) then
            if (ifillctrline==0) write(*,*) "9 Enable filling colors for contour lines"
            if (ifillctrline==1) write(*,*) "9 Set status of filling colors between the contour lines"
        end if
		if (idrawtype==6) then
			if (igrad_arrow==0) write(*,*) "10 Show arrow on the gradient lines"
			if (igrad_arrow==1) write(*,*) "10 Do not show arrow on the gradient lines"
            write(*,*) "11 Set detailed parameters of plotting gradient line"
            if (igrad_arrow==1) then
				if (iinvgradvec==0) write(*,*) "13 Invert gradient vectors"
				if (iinvgradvec==1) write(*,*) "13 Do not invert gradient vectors"
            end if
		else if (idrawtype==7) then
			write(*,"(a,f8.4)") " 10 Set upper limit of absolute value for scaling arrows, current:",cutgradvec
			if (icolorvecfield==0) write(*,*) "11 Map color to arrows"
			if (icolorvecfield==1) write(*,*) "11 Do not map color to arrows"
			if (icolorvecfield==0) write(*,*) "12 Set color for arrow heads" !If color map was set, the color set by user themselves is nulified
			if (iinvgradvec==0) write(*,*) "13 Invert gradient vectors"
			if (iinvgradvec==1) write(*,*) "13 Do not invert gradient vectors"
		end if			
		if (idrawplanevdwctr==0.and.iorbsel2==0.and.allocated(b)) write(*,*) "15 Enable showing the contour line corresponding to vdW surface (rho=0.001)" !meaningless if custom operation is performed
		if (idrawplanevdwctr==1) write(*,*) "15 Disable showing the contour line corresponding to vdW surface (rho=0.001)"
		if (idrawplanevdwctr==1) write(*,*) "16 Set label size, style and color of the contour line of vdW surface"
		if (iatom_on_plane==1) then
			if (idrawintbasple==0) then
				write(*,"(a,f7.3,' Bohr')") " 17 Set distance threshold for showing atom labels, current:",disshowlabel
			else
				write(*,"(a,f6.3,' Bohr')") " 17 Set distance threshold for showing atoms & interbasin paths:",disshowlabel
			end if
		end if
		if (iatmlabtype==1) write(*,*) "18 Change style of atomic labels: Only plot element symbol"
		if (iatmlabtype==2) write(*,*) "18 Change style of atomic labels: Only plot atomic index"
		if (iatmlabtype==3) write(*,*) "18 Change style of atomic labels: Plot both element symbol and atomic index"
        if (iextrema_on_contour==0) then
			write(*,"(a)") " 19 Enable showing extrema of a function on a contour line"
        else
			write(*,"(a)") " 19 Disable showing extrema of a function on a contour line"
        end if
	else if (idrawtype==4.or.idrawtype==5) then !Colored relief map with/without projected color-filled map
		write(*,*) "1 Set color scale range for filling color"
		write(*,"(a,a)") " 2 Toggle drawing mesh on the surface, current: ",drawsurmesh
        write(*,"(a,a)") " 3 Set color transition, current: ",trim(clrtransname(iclrtrans))
	else if (idrawtype==3) then
		continue
	end if
	
	read(*,*) i
	
	!Below are general options for all kinds of plane maps
	if (i==-10) then !Load plane data in another plain text file and operate to current plane data, the plane settings must be identical
		write(*,*) "Input file name, e.g. C:\plane.txt"
		read(*,"(a)") c200tmp
		write(*,*) "How many columns? (4 or 6. The data in the last column will be loaded)"
		read(*,*) ncol
		open(10,file=c200tmp,status="old")
		do i=1,ngridnum1
			do j=1,ngridnum2
				if (ncol==4) then
					read(10,*) tmpv,tmpv,tmpv,planemattmp(i,j)
				else
					read(10,*) tmpv,tmpv,tmpv,tmpv,tmpv,planemattmp(i,j)
				end if
			end do
		end do
		close(10)
		write(*,*) "Which operation? Available operators: +,-,x,/"
		read(*,*) c200tmp(1:1)
		if (c200tmp(1:1)=="+") then
			planemat=planemat+planemattmp
		else if (c200tmp(1:1)=="-") then
			planemat=planemat-planemattmp
		else if (c200tmp(1:1)=="x") then
			planemat=planemat*planemattmp
		else if (c200tmp(1:1)=="/") then
			planemat=planemat/planemattmp
		end if
		write(*,*) "Done!"
	else if (i==-9) then
		if (iatmcontri==0) then
			allocate(planemat_bk(ngridnum1,ngridnum2))
			planemat_bk=planemat
			if (allocated(tmparrint)) deallocate(tmparrint)
			write(*,*) "Input index of the atoms you are interested in, e.g. 2,3,7-10"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,ntmp)
			allocate(tmparrint(ntmp))
			call str2arr(c2000tmp,ntmp,tmparrint)
			write(*,*) "Updating plane data, please wait..."
			do i=1,ngridnum1 !First calculate promolecular density and store it to planemat
				do j=1,ngridnum2
                    call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
					densall=0
					densfrag=0
					do iatm=1,ncenter
						tmpval=calcatmdens(iatm,rnowx,rnowy,rnowz,0)
						densall=densall+tmpval
						if (any(tmparrint==iatm)) densfrag=densfrag+tmpval
					end do
					planemat(i,j)=planemat(i,j)*densfrag/densall
				end do
			end do
			write(*,*) "Done! The data have been updated, you can replot it"
			deallocate(tmparrint)
			iatmcontri=1
		else !Recovery the backed up data
			planemat=planemat_bk
			deallocate(planemat_bk)
			iatmcontri=0
		end if
	else if (i==-8) then
		if (ilenunit2D==1) then
			ilenunit2D=2
		else if (ilenunit2D==2) then
			ilenunit2D=1
		end if
	else if (i==-7) then		
		write(*,*) "Input the value to be multiplied to the current plane data, e.g. 0.3"
		read(*,*) scaleval
		planemat=planemat*scaleval
        gradd1=gradd1*scaleval
        gradd2=gradd2*scaleval
		write(*,*) "Done!"
	else if (i==-6) then
        c200tmp="plane.txt"
        if (iaddprefix==1) call addprefix(c200tmp)
		open(10,file=c200tmp,status="replace")
		do i=1,ngridnum1
			do j=1,ngridnum2
                call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
				if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
					write(10,"(5f10.5,1PE20.10E3)") rnowx*b2a,rnowy*b2a,rnowz*b2a,i*d1*b2a,j*d2*b2a,planemat(i,j)
				else !Plane is vertical, the coordinate in a direction is zero
					write(10,"(3f10.5,1PE20.10E3)") rnowx*b2a,rnowy*b2a,rnowz*b2a,planemat(i,j)
				end if
			end do
		end do
		close(10)
		if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
			write(*,"(a)") " The column 1,2,3 correspond to Cartesian X,Y,Z coordinates, respectively"
			write(*,"(a)") " The column 4,5,6 correspond to X,Y coordinates in the graph and function value, respectively"
		else
			write(*,"(a)") " The column 1,2,3,4 correspond to X,Y,Z and function value, respectively"
		end if
		write(*,"(a)") " Data has been ouputted to "//trim(c200tmp)//" in current folder, length unit is in Angstrom"
		if (idrawtype/=3.and.idrawtype/=4.and.idrawtype/=5) then
			if ( numcp>0.or.(numpath>0.and.imarkpath==1).or.(nple3n1path>0.and.idrawintbasple==1) ) then
				write(*,"(/,a)") " If also outputting critical points and topology/basin paths to plain text files in current folder? (y/n)"
				read(*,*) selectyn
				if (selectyn=='y'.or.selectyn=='Y') then
					iplaneoutall=1 !Global variable, which tells drawplane routine to export topology data
					call drawplane(axlow1,axhigh1,axlow2,axhigh2,drawlowlim,drawuplim)
					iplaneoutall=0
					if (numcp>0) write(*,"(a)") " Critical points have been outputted to planeCP.txt in current folder. The third column is type: 1=(3,-3), 2=(3,-1), 3=(3,+1), 4=(3,+3)"
					if (numpath>0.and.imarkpath==1) write(*,*) "Topology paths have been outputted to planepath.txt in current folder"
					if (nple3n1path>0.and.idrawintbasple==1) write(*,*) "Interbasin paths have been outputted to planeinterbasin.txt in current folder"
					write(*,"(a)") " The first two columns correspond to X,Y coordinates in the graph, the unit is in Angstrom"
				end if
			end if
		end if
	else if (i==-5) then !Return to main menu
		deallocate(planemat,planemattmp)
		if (allocated(planemat_bk)) deallocate(planemat_bk)
		idrawplanevdwctr=0
		iorbsel2=0
        iextrema_on_contour=0
		exit
	else if (i==-4) then
        call saveload2Dplottingsetting(drawlowlim,drawuplim)
	else if (i==-3) then
        if (idrawtype==1.or.(idrawtype==2.and.ifillctrline==1).or.idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
            call plane_othersetting(1)
        else
            call plane_othersetting(0)
        end if
	else if (i==-2) then
		if (idrawtype==1.or.idrawtype==4.or.idrawtype==5) then
			write(*,"(a)") " Input interval between the labels in X, Y and color scale axes, e.g. 1.5,2.0,0.1"
			read(*,*) planestpx,planestpy,planestpz
		else
			write(*,"(a)") " Input interval between the labels in X and Y axes, e.g. 1.5,2.0"
			read(*,*) planestpx,planestpy
		end if
	else if (i==-1) then
		cycle
	else if (i==0) then
		isavepic=1
	end if

	if (idrawtype==1.or.idrawtype==2.or.idrawtype==6.or.idrawtype==7) then
		!Shared options by idrawtype 1,2,6,7
		if (i==8) then
			if (ibond_on_plane==0) then
				ibond_on_plane=1
				write(*,"(a,/)") " Note: The bonding will be empirically determined according to interatomic distance and atomic covalent radii"
				write(*,*) "Use which color for drawing the bonds?"
				call selcolor(iclrindbndlab)
			else
				ibond_on_plane=0
			end if
		else if (i==15) then
			if (idrawplanevdwctr==0) then
				idrawplanevdwctr=1
				write(*,*) "Please wait..."
				!$OMP PARALLEL DO private(ipt,jpt,rnowx,rnowy,rnowz) shared(planemattmp) schedule(dynamic) NUM_THREADS(nthreads)
				do ipt=1,ngridnum1
					do jpt=1,ngridnum2
                        call get2Dgridxyz(ipt,jpt,rnowx,rnowy,rnowz)
						planemattmp(ipt,jpt)=fdens(rnowx,rnowy,rnowz)
					end do
				end do
				!$OMP END PARALLEL DO
				write(*,*) "Done, now you can replot the graph to check effect"
			else if (idrawplanevdwctr==1) then
				idrawplanevdwctr=0
			end if
		else if (i==16) then
			write(*,*) "Input the size of the label, e.g. 30"
			write(*,*) "Note: If you input 0, then the label will not be shown"
			read(*,*) ivdwctrlabsize
			write(*,*) "Select color of the contour line:"
			call selcolor(ivdwclrindctr)
			write(*,*) "Input the width of the contour line, e.g. 10"
			read(*,*) iwidthvdwctr
			write(*,*) "Input length of line segment and interstice"
			write(*,*) "e.g. 1,0 means solid line; 1,10 means DOT; 10,15 means DASH"
			write(*,*) "     10,25 means DASH with larger interstice"
			read(*,*) vdwctrstyle
			write(*,*) "Done, now you can replot the graph to check effect"
		else if (i==17) then
			write(*,*) "Input distance threshold for plotting atomic labels (in Bohr), e.g. 0.5"
            write(*,*) "If you want to input in Angstrom, add ""A"" suffix, e.g. 0.45 A"
			write(*,*) "Note: The default value can be set by ""disshowlabel"" in settings.ini"
            read(*,"(a)") c80tmp
			read(c80tmp,*) disshowlabel
            if (index(c80tmp,'A')/=0.or.index(c80tmp,'a')/=0) disshowlabel=disshowlabel/b2a
			disshowCP=disshowlabel
			disshowpath=disshowlabel
			if (numCP>0.or.numpath>0) write(*,"(a,/)") " Note: The distance threshold for showing CPs/paths has also been changed to this value"
			write(*,"(a)") " If also show labels of the atoms that beyond this criterion as light face type? (y/n)"
			read(*,*) selectyn
			if (selectyn=='y'.or.selectyn=='Y') then
				iatom_on_plane_far=1
			else
				iatom_on_plane_far=0
			end if
		else if (i==18) then
			write(*,*) "1: Only plot element symbol"
			write(*,*) "2: Only plot atomic index"
			write(*,*) "3: Plot both element symbol and atomic index"
			write(*,*) "Note that the default setting can be set by ""iatmlabtype"" in settings.ini"
			read(*,*) iatmlabtype
		else if (i==19) then
			if (iextrema_on_contour==0) then
				write(*,*) "Input isovalue for defining the contour line, e.g. 0.5"
				read(*,*) ctrval_2Dextrema
				write(*,*) "Select the real space function to search for extrema on the contour lines"
				call selfunc_interface(1,ifunc_2Dextrema)
                iextrema_on_contour=1
                if (ifunc_2Dextrema==12) call doinitlibreta(info)
            else
				iextrema_on_contour=0
            end if
		end if
		
		!Options only for idrawtype 1 =====================
		if (idrawtype==1) then 
			if (i==1) then
				write(*,*) "Input lower & upper limit of Z, e.g. -0.3,0.5"
				read(*,*) drawlowlim,drawuplim
			else if (i==2) then
				if (idrawcontour==1) then
					idrawcontour=0
				else if (idrawcontour==0) then
					idrawcontour=1
				end if
			else if (i==3) then  !! Change isovalues of contour line
				call setcontour
			else if (i==4) then
				if (iatom_on_plane==1) then
					iatom_on_plane=0
				else if (iatom_on_plane==0) then
					iatom_on_plane=1
					write(*,*) "Use which color for labelling atoms?"
					call selcolor(iclrindatmlab)
				end if
			else if (i==5) then
				call settopomark
			else if (i==19) then
                call selcolortable
			end if
		!Options only for idrawtype 2,6,7 ========================
		else if (idrawtype==2.or.idrawtype==6.or.idrawtype==7) then
			!General option for idrawtype 2,6,7
			if (i==1) then
				if (iatom_on_plane==1) then
					iatom_on_plane=0
				else if (iatom_on_plane==0) then
					iatom_on_plane=1
					write(*,*) "Use which color?"
					call selcolor(iclrindatmlab)
				end if
			else if (i==2) then
				if (ilabel_on_contour==1) then
					ilabel_on_contour=0
				else if (ilabel_on_contour==0) then
					ilabel_on_contour=1
					write(*,*) "Input label size, e.g. 30"
					read(*,*) ictrlabsize
					write(*,"(a)") " Note: The number of decimal places of labels on contour lines can be set by ""numdigctr"" in settings.ini"
				end if
			else if (i==3) then  !! Change isovalues of contour line
				call setcontour
			else if (i==4) then
				call settopomark
			else if (i==5) then
				if (idrawcontour==1) then
					idrawcontour=0
				else if (idrawcontour==0) then
					idrawcontour=1
				end if
			else if (i==6) then !Generate paths from (3,-1)
				if (idrawintbasple==0) then	
					do icp=1,numcp
						if (CPtype(icp)==2) then
							cpx=CPpos(1,icp)
							cpy=CPpos(2,icp)
							cpz=CPpos(3,icp)
							if (plesel==1) then
								if (abs(cpz-orgz2D) > disshowlabel) cycle
							else if (plesel==2) then
								if (abs(cpy-orgy2D) > disshowlabel) cycle
							else if (plesel==3) then
								if (abs(cpx-orgx2D) > disshowlabel) cycle
							else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
								call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,cpx,cpy,cpz,prjx,prjy,prjz)
								if ( (cpx-prjx)**2+(cpy-prjy)**2+(cpz-prjz)**2 > disshowlabel**2) cycle
							end if
							nple3n1path=nple3n1path+1
							cp2ple3n1path(icp)=nple3n1path !Default is zero, means this CP is not on the given plane and has no corresponding interbasin path
						end if
					end do
					if (nple3n1path>0) then
						idrawintbasple=1
						write(*,"(' Found',i8,' (3,-1) CPs in the plane')") nple3n1path
						allocate(ple3n1path(3,n3n1plept,2,nple3n1path))
						write(*,*) "Generating interbasin paths from (3,-1) CPs, Please wait..."
						!$OMP PARALLEL DO SHARED(numcp) PRIVATE(icp) schedule(dynamic) NUM_THREADS(nthreads)
						do icp=1,numcp
							if (cp2ple3n1path(icp)/=0) call gen3n1plepath(ifunctopo,icp,cp2ple3n1path(icp))
 							! write(*,"('Finished the in-plane path generation from (3,-1)',i8)") icp
						end do
						!$OMP END PARALLEL DO
					else
						write(*,*) "No (3,-1) CP is closed to the plane"
					end if
				else if (idrawintbasple==1) then
					idrawintbasple=0
					nple3n1path=0
					cp2ple3n1path=0
					if (allocated(ple3n1path)) deallocate(ple3n1path)
				end if
			else if (i==7) then
				write(*,*) "Input stepsize (Bohr) and maximal iterations, e.g. 0.01,200"
				write(*,"(a,f8.5,',',i6)") " Current values:",ple3n1pathstpsiz,n3n1plept
				read(*,*) ple3n1pathstpsiz,n3n1plept
            else if (i==9) then
                if (ifillctrline==0) then
                    ifillctrline=1
                else
                    do while(.true.)
                        write(*,*)
                        write(*,*) "0 Return"
                        write(*,*) "1 Disable coloring contour lines"
                        write(*,"(a,1PE13.5,a,1PE13.5)") " 2 Set lower and upper limits of filling, current:",drawlowlim," to",drawuplim
                        write(*,"(a,a)") " 3 Set color transition, current: ",trim(clrtransname(iclrtrans))
                        if (ishowclrfill_bar==1) write(*,*) "4 Toggle showing color bar, current: Yes"
                        if (ishowclrfill_bar==0) write(*,*) "4 Toggle showing color bar, current: No"
                        if (ishowclrfill_bar==1) write(*,*) "5 Set label interval of color bar"
                        read(*,*) itmp
                        if (itmp==0) then
                            exit
                        else if (itmp==1) then
                            ifillctrline=0
                        else if (itmp==2) then
                            write(*,*) "Input lower and upper limits of color scale, e.g. -0.3,0.3"
                            read(*,*) drawlowlim,drawuplim
                        else if (itmp==3) then
                            call selcolortable
                        else if (itmp==4) then
                            if (ishowclrfill_bar==1) then
                                ishowclrfill_bar=0
                            else
                                ishowclrfill_bar=1
                            end if
                        else if (itmp==5) then
                            write(*,*) "Input label interval, e.g. 0.2"
                            read(*,*) planestpz
                        end if
                    end do
                end if
			end if
			
			!Option only for idrawtype 6 ===========
			if (idrawtype==6) then
            
				if (i==10) then
					if (igrad_arrow==1) then
						igrad_arrow=0
					else if (igrad_arrow==0) then
						igrad_arrow=1
					end if
				else if (i==11) then
                    do while(.true.)
                        write(*,*)
                        write(*,*) "   ------------- Set parameters for plotting gradient lines -------------"
                        write(*,*) "0 Return"
			            write(*,"(a,f8.4)") " 1 Set integration step of gradient lines, current:",gradplotstep
			            write(*,"(a,f8.4)") " 2 Set interstice between gradient lines, current:",gradplotdis
			            write(*,"(a,f8.4)") " 3 Set test value of drawing a new gradient line, current:",gradplottest
                        write(*,*) "4 Set integration method, current: "//stream_intmethod
			            write(*,*) "5 Set color width of gradient lines, current: "//colorname(iclrindgradline)
			            write(*,"(a,i3)") " 6 Set line width of gradient lines, current:",iwidthgradline
                        read(*,*) itmp
                        if (itmp==0) then
                            exit
                        else if (itmp==1) then
					        write(*,*) "Input a value, the default value is 0.002"
					        write(*,"(a)") " Note: The smaller the value, the smoother the gradient lines, but higher the cost"
					        read(*,*) gradplotstep
				        else if (itmp==2) then
					        write(*,*) "Input a value, the default value is 0.01"
					        write(*,*) "Note: The larger the value, the sparser the gradient lines"
					        read(*,*) gradplotdis
				        else if (itmp==3) then
					        write(*,*) "Input a value, default value is 0.2"
					        write(*,"(a)") " Note: The smaller the value, the more gradient lines will be emitted for maxima"
					        read(*,*) gradplottest
				        else if (itmp==4) then
                            write(*,*) "1 Euler"
                            write(*,*) "2 RK2 (default)"
                            write(*,*) "3 RK4 (most robust)"
                            read(*,*) iint
                            if (iint==1) stream_intmethod="EULER"
                            if (iint==2) stream_intmethod="RK2"
                            if (iint==3) stream_intmethod="RK4"
				        else if (itmp==5) then
					        write(*,*) "Use which color?"
					        call selcolor(iclrindgradline)
				        else if (itmp==6) then
					        write(*,*) "Input line width, e.g. 5  (default value is 1)"
					        read(*,*) iwidthgradline
                        end if
                    end do
				else if (i==13) then
					if (iinvgradvec==1) then
						iinvgradvec=0
					else if (iinvgradvec==0) then
						iinvgradvec=1
					end if
				end if
			!Option only for idrawtype 7 ============
			else if (idrawtype==7) then
				if (i==10) then
					write(*,*) "Input a value, e.g. 0.3"
					read(*,*) cutgradvec
				else if (i==11) then
					if (icolorvecfield==1) then
						icolorvecfield=0
					else if (icolorvecfield==0) then
						icolorvecfield=1
					end if
				else if (i==12) then
					write(*,*) "Input color index, e.g. 1=black, 50=blue, 150=green, 250=red"
					read(*,*) vecclrind
				else if (i==13) then
					if (iinvgradvec==1) then
						iinvgradvec=0
					else if (iinvgradvec==0) then
						iinvgradvec=1
					end if
				end if
			end if
		end if
		
	!Options for idrawtype 4,5 =================
	else if (idrawtype==4.or.idrawtype==5) then
		if (i==1) then
			write(*,*) "Input lower & upper limits of color scale for shading surface, e.g. -0.1,0.3"
			write(*,"(a,f14.6,a,f14.6)") " Present value: from",surcolorzmin," to",surcolorzmax
			read(*,*) surcolorzmin,surcolorzmax
			if (idrawtype==5) then
				write(*,*) "Input lower & upper limits of color scale for projected map, e.g. -0.1,0.3"
				write(*,"(a,f14.6,a,f14.6)") " Present value: from",drawlowlim," to",drawuplim
				read(*,*) drawlowlim,drawuplim
			end if
		else if (i==2) then
			if (drawsurmesh=="ON ") then
				drawsurmesh="OFF"
			else
				drawsurmesh="ON"
			end if
		else if (i==3) then
            CALL selcolortable
		end if
	else if (idrawtype==3) then !No options for map type 3 currently
		continue
	end if
end	do

end subroutine




!!---------- Define a plane above/below the plane consisting of specific atoms    
subroutine define_parallel_plane
use defvar
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c2000tmp*2000
integer,allocatable :: tmparr(:)
real*8 pt1(3),pt2(3),pt3(3),vec1(3),vec2(3),vectmp(3),cenpos(3)

write(*,*) "Input indices of the atoms to define the fitting plane and geometric center"
write(*,*) "e.g. 3,8,10-15,19"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,ntmp)
allocate(tmparr(ntmp))
call str2arr(c2000tmp,ntmp,tmparr)
call ptsfitplane(tmparr,ntmp,planeA,planeB,planeC,rnouse,rmsfit)
rnorm=dsqrt(planeA**2+planeB**2+planeC**2)
write(*,"(' The unit normal vector is',3f14.8)") planeA/rnorm,planeB/rnorm,planeC/rnorm
write(*,*)
write(*,*) "Input distance between plotting plane and geometric center in Angstrom, e.g. 1"
write(*,*) "Positive and negative values corresponding to above and below the plane"
read(*,*) pledist
pledist=pledist/b2a

write(*,*) "Input length of the plotting plane in Angstrom, e.g. 5.5"
write(*,"(a)") " Hint: Usually 6 Angstrom is adequate for studying a six-membered ring. For a larger region, use a larger value"
read(*,*) plelen
plelen=plelen/b2a

!Project arbitrary three points onto the Ax+By+Cz=0 plane
call pointprjpleABCD(planeA,planeB,planeC,0D0,3.5D0,1.2D0,-5.2D0,pt1(1),pt1(2),pt1(3))
call pointprjpleABCD(planeA,planeB,planeC,0D0,2.5D0,0.7D0,-3.3D0,pt2(1),pt2(2),pt2(3))
call pointprjpleABCD(planeA,planeB,planeC,0D0,6.5D0,-1.2D0,2.7D0,pt3(1),pt3(2),pt3(3))

!Make vec2 (3-1) ortho to vec1 (2-1) using Schmit method, namely
!3
!|
!1---2
vec1(:)=pt2(:)-pt1(:)
vec2(:)=pt3(:)-pt1(:)
schmit=dot_product(vec1,vec2)/sum(vec1**2)
vec2(:)=vec2(:)-schmit*vec1(:)

!Normalization to expected length
vec1(:)=vec1(:)/dsqrt(sum(vec1**2))*plelen
vec2(:)=vec2(:)/dsqrt(sum(vec2**2))*plelen

!Determine the position above ring center
cenpos(1)=sum(a(tmparr(:))%x)/ntmp
cenpos(2)=sum(a(tmparr(:))%y)/ntmp
cenpos(3)=sum(a(tmparr(:))%z)/ntmp
vectmp(1)=planeA
vectmp(2)=planeB
vectmp(3)=planeC
cenpos(:)=cenpos(:)+pledist*vectmp/dsqrt(sum(vectmp**2))

!Determine position of three points
vectmp(:)=(-vec1(:)-vec2(:))/2D0
a1x=cenpos(1)+vectmp(1)
a1y=cenpos(2)+vectmp(2)
a1z=cenpos(3)+vectmp(3)
a2x=a1x+vec1(1)
a2y=a1y+vec1(2)
a2z=a1z+vec1(3)
a3x=a1x+vec2(1)
a3y=a1y+vec2(2)
a3z=a1z+vec2(3)
!write(*,"(' Point 1:'3f12.6,' Angstrom')") a1x*b2a,a1y*b2a,a1z*b2a
!write(*,"(' Point 2:'3f12.6,' Angstrom')") a2x*b2a,a2y*b2a,a2z*b2a
!write(*,"(' Point 3:'3f12.6,' Angstrom')") a3x*b2a,a3y*b2a,a3z*b2a
a4x=cenpos(1)-vectmp(1)
a4y=cenpos(2)-vectmp(2)
a4z=cenpos(3)-vectmp(3)
write(*,*) "Command of drawing the plotting plane in VMD:"
write(*,"('draw triangle {',3f8.3,'} {',3f8.3,'} {',3f8.3,'}')") a1x*b2a,a1y*b2a,a1z*b2a,a2x*b2a,a2y*b2a,a2z*b2a,a3x*b2a,a3y*b2a,a3z*b2a
write(*,"('draw triangle {',3f8.3,'} {',3f8.3,'} {',3f8.3,'}')") a2x*b2a,a2y*b2a,a2z*b2a,a3x*b2a,a3y*b2a,a3z*b2a,a4x*b2a,a4y*b2a,a4z*b2a
write(*,"('draw material Transparent')")
end subroutine




!!-------- Interface of changing other settings
!If this plane plotting involve Z-axis, iZaxis should be 1, otherwise 0
subroutine plane_othersetting(iZaxis)
use defvar
implicit real*8 (a-h,o-z)
integer iZaxis
do while(.true.)
    write(*,*)
    write(*,*) "               ------------- Other plotting settings -------------"
    write(*,*) "0 Return"
    if (itickreverse==0) write(*,*) "1 Toggle reversing ticks, current: No"
    if (itickreverse==1) write(*,*) "1 Toggle reversing ticks, current: Yes"
    write(*,*) "2 Set number of decimal places of tick labels"
    write(*,"(a,2i6)") " 3 Set pixel for width and height of the exported figure, current:",graph2Dwidth,graph2Dheight
	write(*,"(a,i3)") " 4 Set the number of ticks between the labels, current:",iticks-1
    write(*,"(a,i3)") " 5 Set text size of tick labels, current:",plane_axistextsize
    write(*,"(a,i3)") " 6 Set text size of axis names, current:",plane_axisnamesize
	if (iatmlabtype==1) write(*,*) "7 Change style of atomic labels: Only plot element symbol"
	if (iatmlabtype==2) write(*,*) "7 Change style of atomic labels: Only plot atomic index"
	if (iatmlabtype==3) write(*,*) "7 Change style of atomic labels: Plot both element symbol and atomic index"
    write(*,"(a,i4)") " 8 Set size of atomic labels, current:",pleatmlabsize
    write(*,"(a,i3)") " 9 Set thick of bonds, current:",bondthick2D
	write(*,*) "10 Set format of exporting image file, current: "//graphformat
	if (idrawtype==2) write(*,"(a,i3)") " 11 Set number of decimal places of isovalue labels on contours, current:",numdigctr
    if (idrawtype==1.and.idrawcontour==1) then
    		if (ilabel_on_contour==0) write(*,*) "12 Enable showing isovalue on contour lines"
		if (ilabel_on_contour==1) write(*,*) "12 Disable showing isovalue on contour lines"
    end if
    read(*,*) isel2
    
    if (isel2==0) then
        exit
    else if (isel2==1) then
		if (itickreverse==0) then
			itickreverse=1
		else
			itickreverse=0
		end if
    else if (isel2==2) then
        write(*,*) "Note: The default value can be set by ""numdigxyz"" parameter in settings.ini"
        write(*,*)
        write(*,*) "Input number of digits after the decimal point for X axis, e.g. 2"
        write(*,"(' Current value:',i3)") numdigx
        read(*,*) numdigx
        write(*,*) "Input number of digits after the decimal point for Y axis, e.g. 2"
        write(*,"(' Current value:',i3)") numdigy
        read(*,*) numdigy
        if (iZaxis==1) then
            write(*,*) "Input number of digits after the decimal point for Z axis, e.g. 3"
            write(*,"(' Current value:',i3)") numdigz
            read(*,*) numdigz
        end if
        write(*,*) "Done!"
    else if (isel2==3) then
        write(*,*) "Note: The default value can be set by ""graph2Dsize"" parameter in settings.ini"
        write(*,*)
        write(*,*) "Input pixel for width, e.g. 1500"
        write(*,"(' Current value:',i5)") graph2Dwidth
        read(*,*) graph2Dwidth
        write(*,*) "Input pixel for height, e.g. 1200"
        write(*,"(' Current value:',i5)") graph2Dheight
        write(*,*) "Note: If you input 0, then height will be set to",int(dfloat(graph2Dwidth)/5*4)
        read(*,*) graph2Dheight
        if (graph2Dheight==0) graph2Dheight=dfloat(graph2Dwidth)/5*4
        write(*,*) "Done!"
    else if (isel2==4) then
		write(*,*) "How many ticks between labels do you want? e.g. 2"
		read(*,*) iticks
		iticks=iticks+1
    else if (isel2==5) then
        write(*,*) "Input text size of ticks on the axes, e.g. 45"
        write(*,"(' Current value:',i3)") plane_axistextsize
        read(*,*) plane_axistextsize
    else if (isel2==6) then
        write(*,*) "Input text size of axis names, e.g. 40"
        write(*,"(' Current value:',i3)") plane_axisnamesize
        read(*,*) plane_axisnamesize
    else if (isel2==7) then
		write(*,*) "1: Only plot element symbol"
		write(*,*) "2: Only plot atomic index"
		write(*,*) "3: Plot both element symbol and atomic index"
		write(*,*) "Note: The default value can be set by ""iatmlabtype"" in settings.ini"
		read(*,*) iatmlabtype
    else if (isel2==8) then
        write(*,*) "Input size of atomic labels, e.g. 40"
        write(*,"(' Current value:',i3)") pleatmlabsize
		write(*,*) "Note: The default value can be set by ""pleatmlabsize"" in settings.ini"
        read(*,*) pleatmlabsize
    else if (isel2==9) then
        write(*,*) "Input thick of bonds, e.g. 20"
        write(*,"(' Current value:',i3)") bondthick2D
        read(*,*) bondthick2D
    else if (isel2==10) then
        call setgraphformat
    else if (isel2==11) then
        write(*,*) "Input number of decimal places on labels, e.g. 4"
        write(*,"(' Current value:',i3)") numdigctr
        write(*,*) "Note: The default value can be set by ""numdigctr"" parameter in settings.ini"
        read(*,*) numdigctr
    else if (isel2==12) then
		if (ilabel_on_contour==0) then
			ilabel_on_contour=1
			write(*,*) "Input label size, e.g. 30"
			read(*,*) ictrlabsize
			write(*,"(a)") " Note: The number of decimal places of labels on contour lines can be set by ""numdigctr"" in settings.ini"
        else
			ilabel_on_contour=0
        end if
    end if
end do
end subroutine






!===============================================================================================
!===============================================================================================
!---------- Study various real space functions in a spatial region (i.e. 3 dimension) ----------
!===============================================================================================
!===============================================================================================
subroutine study3dim
use defvar
use util
use functions
use GUI
implicit real*8 (a-h,o-z)
character c200tmp*200,c2000tmp*2000,outcubfile*200
real*8 :: tmpvec(3)
integer,allocatable :: tmparrint(:)
real*8,allocatable :: cubmat_bk(:,:,:)
!Wsed when generate grid data based on external grid data by interpolation
real*8 gridv1_ext(3),gridv2_ext(3),gridv3_ext(3),orgx_ext,orgy_ext,orgz_ext
integer nx_ext,ny_ext,nz_ext
real*8 gridv1_new(3),gridv2_new(3),gridv3_new(3),orgx_new,orgy_new,orgz_new
integer nx_new,ny_new,nz_new

ncustommap=0 !Clean custom operation setting that possibly defined by other modules
if (allocated(custommapname)) deallocate(custommapname)
if (allocated(customop)) deallocate(customop)
ipromol=0
	
do while(.true.)
    write(*,*) "-10 Return to main menu"
	write(*,*) "-2 Obtain deformation property"
	write(*,*) "-1 Obtain promolecule property"
	write(*,*) "0 Set custom operation"
	call selfunc_interface(1,ifuncsel)
	
	if (ifuncsel==-10) then
		return
	else if (ifuncsel==-2) then
		call setPromol
		customop='-'
	else if (ifuncsel==-1) then
		ipromol=1 !Special case, obtain promolecular property
		call setPromol
		customop='+'
	else if (ifuncsel==0) then
		call customplotsetup
	else if (ifuncsel==111) then !Calculate Becke weighting function
		write(*,*) "Input indices of two atoms to calculate Becke overlap weight, e.g. 1,4"
		write(*,*) "or input index of an atom and zero to calculate Becke atomic weight, e.g. 5,0"
		read(*,*) iatmbecke1,iatmbecke2
		exit
	else if (ifuncsel==112) then !Calculate Hirshfeld weighting function
		write(*,*) "Input index of the atoms that you want to calculate Hirshfeld weight"
		write(*,*) "e.g. 2,3,7-10"
		read(*,"(a)") c2000tmp
		call str2arr(c2000tmp,ntmp)
		allocate(tmparrint(ntmp))
		call str2arr(c2000tmp,ntmp,tmparrint)
		write(*,"(a)") " How to generate the atomic densities that used in the calculation of Hirshfeld weight?"
		write(*,*) "1 Based on atomic .wfn files"
		write(*,*) "2 Based on built-in atomic densities (see Appendix 3 of the manual for detail)"
		read(*,*) iHirshdenstype
		if (iHirshdenstype==1) call setpromol
		exit
	else
		exit
	end if
end do

if (ifuncsel==100.and.(iuserfunc==-1.or.iuserfunc==-3)) then !Special case, interpolation based on grid data, backup grid setting of external grid
	gridv1_ext(:)=gridv1(:)
    gridv2_ext(:)=gridv2(:)
    gridv3_ext(:)=gridv3(:)
    orgx_ext=orgx;orgy_ext=orgy;orgz_ext=orgz
    nx_ext=nx;ny_ext=ny;nz_ext=nz
end if

call setgrid(1,igridsel)

if (igridsel==100) then !Calculate value on a set of points loaded from external file, and then write into the last column of the file
	if (allocated(extpttmp)) deallocate(extpttmp)
	if (ncustommap/=0) allocate(extpttmp(numextpt)) !!! temp file for difference cube
	if (ifuncsel/=4) call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
	icustom=0
	extpt(:,4)=0D0
	
	call walltime(iwalltime1)
	
	if (ipromol==1) goto 501 !Calculate promolecular property, so skip the first time calculation (namely for the whole system)
500 continue
	!$OMP PARALLEL DO SHARED(extpt) PRIVATE(iextpt) schedule(dynamic) NUM_THREADS(nthreads)
	do iextpt=1,numextpt !Calculate function value
		extpt(iextpt,4)=calcfuncall(ifuncsel,extpt(iextpt,1),extpt(iextpt,2),extpt(iextpt,3))
	end do
	!$OMP END PARALLEL DO
501	if (ncustommap/=0) then !Cycling will stop when all the file have been dealed
		if (icustom==0) then
	!Note: For promolecular property, x,y,z hasn't been saved in extpt at first time, while after calculation of atoms, extpt already has %x,%y,%z
			extpttmp(:)=extpt(:,4) !first time
		else if (icustom/=0) then !not first time
			if (customop(icustom)=='+') extpttmp(:)=extpttmp(:)+extpt(:,4)
			if (customop(icustom)=='-') extpttmp(:)=extpttmp(:)-extpt(:,4)
			if (customop(icustom)=='x'.or.customop(icustom)=='*') extpttmp(:)=extpttmp(:)*extpt(:,4)
			if (customop(icustom)=='/') extpttmp(:)=extpttmp(:)/extpt(:,4)
		end if
		if (icustom/=ncustommap) then
			icustom=icustom+1
			filename=custommapname(icustom)
            call savePBCinfo
			call dealloall(0)
			write(*,"(' Loading:  ',a)") trim(filename)
			call readinfile(filename,1)
            call loadPBCinfo
			if (ifuncsel/=4) call delvirorb(0)
			!Input the MO index for current file. Since the MO index may be not the same as the first loaded one
			if (ifuncsel==4) then
				write(*,"(' Input index of the orbital to be calculated for ',a,', e.g. 3')") trim(filename)
				read(*,*) iorbsel
			end if
			goto 500
		else if (icustom==ncustommap) then !last time
			extpt(:,4)=extpttmp(:)
			call dealloall(0)
			write(*,"(' Reloading:  ',a)") trim(firstfilename)
			call readinfile(firstfilename,1)
		end if
	end if
    
	call delvirorb_back(1) !delvirorb may have taken effect, now restore to previous wavefunction
	call walltime(iwalltime2)
	write(*,"(' Calculation is finished, took up wall clock time',i10,' s')") iwalltime2-iwalltime1
	
	write(*,"(a)") " Output the points with function values to which file? e.g. C:\ltwd.txt"
	read(*,"(a)") c200tmp
	open(10,file=c200tmp,status="replace")
	write(10,"(i10)") numextpt
	do iextpt=1,numextpt
		write(10,"(3f13.7,E20.10)") extpt(iextpt,1:3),extpt(iextpt,4)
	end do
	close(10)
	write(*,"(a)") " Done! In this file the first line is the number of points, Column 1~4 correspond to X,Y,Z coordinates and function values, respectively. All units are in a.u."

else !Calculate grid data
	if (ifuncsel==100.and.(iuserfunc==-1.or.iuserfunc==-3)) then !Special case, interpolation based on grid data. To avoid overwrite cubmat during interpolation, we use special treatment
		!Backup new grid setting
		gridv1_new(:)=gridv1(:);gridv2_new(:)=gridv2(:);gridv3_new(:)=gridv3(:)
		orgx_new=orgx;orgy_new=orgy;orgz_new=orgz
		nx_new=nx;ny_new=ny;nz_new=nz
        !Restore setting of external grid data, so we can normally use interpolation function
		gridv1(:)=gridv1_ext(:);gridv2(:)=gridv2_ext(:);gridv3(:)=gridv3_ext(:)
        orgx=orgx_ext;orgy=orgy_ext;orgz=orgz_ext
        nx=nx_ext;ny=ny_ext;nz=nz_ext
        if (allocated(cubmattmp)) deallocate(cubmattmp)
        allocate(cubmattmp(nx_new,ny_new,nz_new))
        ifinish=0;ishowprog=1
		ntmp=floor(ny_new*nz_new/100D0)
		!$OMP PARALLEL DO SHARED(cubmattmp,ifinish,ishowprog) PRIVATE(i,j,k,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
        do k=1,nz_new
			do j=1,ny_new
				do i=1,nx_new
					tmpx = orgx_new + gridv1_new(1)*(i-1) + gridv2_new(1)*(j-1) + gridv3_new(1)*(k-1)
					tmpy = orgy_new + gridv1_new(2)*(i-1) + gridv2_new(2)*(j-1) + gridv3_new(2)*(k-1)
					tmpz = orgz_new + gridv1_new(3)*(i-1) + gridv2_new(3)*(j-1) + gridv3_new(3)*(k-1)
					cubmattmp(i,j,k)=calcfuncall(ifuncsel,tmpx,tmpy,tmpz)
				end do
				if (ntmp/=0) then
					!$OMP CRITICAL
					ifinish=ifinish+1
					ishowprog=mod(ifinish,ntmp)
					if (ishowprog==0) call showprog(floor(100D0*ifinish/(ny_new*nz_new)),100)
					!$OMP END CRITICAL
                end if
			end do
		end do
		!$OMP END PARALLEL DO
		if (ishowprog/=0) call showprog(100,100)
        !Restore setting of new grid data
		gridv1(:)=gridv1_new(:);gridv2(:)=gridv2_new(:);gridv3(:)=gridv3_new(:)
		orgx=orgx_new;orgy=orgy_new;orgz=orgz_new
		nx=nx_new;ny=ny_new;nz=nz_new
        deallocate(cubmat)
        allocate(cubmat(nx,ny,nz))
        cubmat=cubmattmp
        deallocate(cubmattmp)
    else !Common case
		if (allocated(cubmat)) deallocate(cubmat)
		if (allocated(cubmattmp)) deallocate(cubmattmp)
		allocate(cubmat(nx,ny,nz))
		if (ncustommap/=0) allocate(cubmattmp(nx,ny,nz)) !!! temp file for difference cube
		if (ifuncsel/=4) call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
    end if
	
    if (ifuncsel==100.and.(iuserfunc==-1.or.iuserfunc==-3)) then
		continue
	else if (ifuncsel==111) then !Becke's weight
		do k=1,nz
			do j=1,ny
				do i=1,nx
                    call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					cubmat(i,j,k)=beckewei(tmpx,tmpy,tmpz,iatmbecke1,iatmbecke2)
				end do
			end do
		end do
	else if (ifuncsel==112) then !Hirshfeld weight
		call genhirshcubewei(tmparrint,size(tmparrint),iHirshdenstype)
		ncustommap=0
	else if (ifuncsel==120) then !Calculate and output three components of Steric force to plain text file
		open(20,file="stericforce.txt",status="replace")
		do k=1,nz
			do j=1,ny
				do i=1,nx
                    call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					call stericderv(tmpx,tmpy,tmpz,tmpvec)
					write(20,"(7f12.6)") tmpx*b2a,tmpy*b2a,tmpz*b2a,-tmpvec,dsqrt(sum(tmpvec**2))
				end do
			end do
			call showprog(k,nz)
		end do
		close(20)
		write(*,*) "Done, the results have been outputted to stericforce.txt in current folder"
		write(*,"(a)") " Columns 1,2,3 correspond to X,Y,Z coordinates, 4,5,6 correspond to &
        &steric force component in X,Y,Z. The last column denotes magnitude of steric force"
		write(*,*)
		read(*,*)
	else if (ifuncsel==112) then !Hirshfeld weight
		call genhirshcubewei(tmparrint,size(tmparrint),iHirshdenstype)
		ncustommap=0
    else if (ifuncsel==100.and.(iuserfunc==57.or.iuserfunc==58.or.iuserfunc==59)) then !Calculate g1,g2,g3 terms defined by Shubin, they rely on rho_0
        call g1g2g3grid
		ncustommap=0
    else if (ifuncsel==504) then !2nd relative Onicescu information
		call genentrocub(4)
		ncustommap=0
    else if (ifuncsel==505) then !3rd relative Onicescu information
		call genentrocub(5)
		ncustommap=0
	else !Common case
		icustom=0
		if (ipromol==1) then
			cubmat=0D0
			goto 511 !Calculate promolecular property, so skip the first time calculation (namely for the whole system)
        end if
	510	call savecubmat(ifuncsel,0,iorbsel) !Calculate data to cubmat matrix
	511	if (ncustommap/=0) then !Calculate data for custom cube, cycling stop when all the file have been dealed
			if (icustom==0) then !first time
			!Note: For promolecular property, x,y,z has not been saved in cubmat at first time, while after calculation of atoms, cubmat already has %x,%y,%z
				cubmattmp=cubmat
			else !Not first time
				if (customop(icustom)=='+') cubmattmp=cubmattmp+cubmat
				if (customop(icustom)=='-') cubmattmp=cubmattmp-cubmat
				if (customop(icustom)=='x'.or.customop(icustom)=='*') cubmattmp=cubmattmp*cubmat
				if (customop(icustom)=='/') cubmattmp=cubmattmp/cubmat
			end if
			if (icustom/=ncustommap) then !Not last time
				icustom=icustom+1
				filename=custommapname(icustom)
                call savePBCinfo !In order to keep PBC status of atom identical to the whole system, saving PBC of the whole system, and then apply to atom after loading atomic wfn file
				call dealloall(0)
				write(*,"(' Loading:  ',a)") trim(filename)
				call readinfile(filename,1)
                call loadPBCinfo
				if (ifuncsel/=4) call delvirorb(0)
				!Input the MO index for current file. Since the MO index may be not the same as the first loaded one
				if (ifuncsel==4) then
					write(*,"(' Input index of the orbital to be calculated for ',a,', e.g. 3')") trim(filename)
					read(*,*) iorbsel
				end if
				goto 510
			else if (icustom==ncustommap) then !last time
				cubmat=cubmattmp
				call dealloall(0)
				write(*,"(' Reloading:  ',a)") trim(firstfilename)
				call readinfile(firstfilename,1)
			end if
		end if
	end if
    call delvirorb_back(1) !delvirorb may have taken effect, now restore to previous wavefunction
    
    call calc_dvol(dvol)
	outcubfile="griddata.cub" !General name
	if (ifuncsel==1) then
		outcubfile="density.cub"
		dipx=0;dipy=0;dipz=0
		do k=1,nz
			do j=1,ny
				do i=1,nx
                    call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					dipx=dipx-cubmat(i,j,k)*tmpx
					dipy=dipy-cubmat(i,j,k)*tmpy
					dipz=dipz-cubmat(i,j,k)*tmpz
				end do
			end do
		end do
        if (nEDFprims/=0) then
		    dipx=sum(a%index*a%x)+dipx*dvol
		    dipy=sum(a%index*a%y)+dipy*dvol
		    dipz=sum(a%index*a%z)+dipz*dvol
        else
		    dipx=sum(a%charge*a%x)+dipx*dvol
		    dipy=sum(a%charge*a%y)+dipy*dvol
		    dipz=sum(a%charge*a%z)+dipz*dvol
        end if
        if (ifPBC==0) then
			write(*,*)
			write(*,*) "Electric dipole moment estimated by integrating electron density"
			write(*,"(' X component:    ',f14.6,' a.u.',f14.6,' Debye')") dipx,dipx*au2debye
			write(*,"(' Y component:    ',f14.6,' a.u.',f14.6,' Debye')") dipy,dipy*au2debye
			write(*,"(' Z component:    ',f14.6,' a.u.',f14.6,' Debye')") dipz,dipz*au2debye
			write(*,"(' Total magnitude:',f14.6,' a.u.',f14.6,' Debye')") dsqrt(dipx**2+dipy**2+dipz**2),dsqrt(dipx**2+dipy**2+dipz**2)*au2debye
        end if
	else if (ifuncsel==2) then
		outcubfile="gradient.cub"
	else if (ifuncsel==3) then
		outcubfile="laplacian.cub"
	else if (ifuncsel==4) then
		outcubfile="MOvalue.cub"
	else if (ifuncsel==5) then
		outcubfile="spindensity.cub"
		sur_value=0.02D0
	else if (ifuncsel==6) then
		outcubfile="K(r).cub"
	else if (ifuncsel==7) then
		outcubfile="G(r).cub"
	else if (ifuncsel==8) then
		outcubfile="nucleiesp.cub"
	else if (ifuncsel==9) then
		outcubfile="ELF.cub"
		sur_value=0.7D0
	else if (ifuncsel==10) then
		outcubfile="LOL.cub"
		sur_value=0.5D0
	else if (ifuncsel==11) then
		outcubfile="infoentro.cub"
	else if (ifuncsel==12) then
		outcubfile="totesp.cub"
	else if (ifuncsel==13) then
		sur_value=0.5D0
		outcubfile="RDG.cub"
	else if (ifuncsel==14) then
		sur_value=0.4D0
		outcubfile="RDGprodens.cub"
	else if (ifuncsel==15) then
		outcubfile="signlambda2rho.cub"
	else if (ifuncsel==16) then
		outcubfile="signlambda2rhoprodens.cub"
	else if (ifuncsel==17) then
		outcubfile="fermihole.cub"
	else if (ifuncsel==18) then
		outcubfile="avglocion.cub"
	else if (ifuncsel==19) then
		outcubfile="srcfunc.cub"
	else if (ifuncsel==20) then
		outcubfile="EDR.cub"
	else if (ifuncsel==21) then
		outcubfile="EDRDmax.cub"
	else if (ifuncsel==22) then
		outcubfile="Delta_g.cub"
	else if (ifuncsel==24) then
		sur_value=1D0
		outcubfile="IRI.cub"
	else if (ifuncsel==25) then
		sur_value=1D0
		outcubfile="vdWpot.cub"
	else if (ifuncsel==44) then
		sur_value=0.005D0
		outcubfile="orbdens.cub"
	else if (ifuncsel==100) then
		outcubfile="userfunc.cub"
	else if (ifuncsel==111) then
		outcubfile="Becke.cub"
	else if (ifuncsel==112) then
		outcubfile="Hirshfeld.cub"
	end if

	temp=minval(cubmat)
	call findvalincub(cubmat,temp,i,j,k)
    call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
	write(*,"(/,' The minimum is',E16.8,' at',3f10.5,' Bohr')") temp,tmpx,tmpy,tmpz
	temp=maxval(cubmat)
	call findvalincub(cubmat,temp,i,j,k)
    call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
	write(*,"(' The maximum is',E16.8,' at',3f10.5,' Bohr')") temp,tmpx,tmpy,tmpz
	write(*,"(' Summing up all value and multiply differential element:')") 
	write(*,*) sum(cubmat)*dvol
	write(*,"(' Summing up positive value and multiply differential element:')")
	write(*,*) sum(cubmat,mask=cubmat>0)*dvol
	write(*,"(' Summing up negative value and multiply differential element:')")
	write(*,*) sum(cubmat,mask=cubmat<0)*dvol

	!Reinitialize plot parameter
	bondcrit=1.15D0
	textheigh=38D0
	ratioatmsphere=1D0
	bondradius=0.2D0
	ishowatmlab=1
	ishowaxis=1
	idrawmol=1
	isosurshowboth=1
	ishowdatarange=0
	
	do while(.true.)
		write(*,*)
        call menutitle("Post-processing menu",10,1)
		write(*,*) "-1 Show isosurface graph"
		write(*,*) "0 Return to main menu"
		write(*,*) "1 Save graph of isosurface to file in current folder"
		write(*,*) "2 Export data to a Gaussian-type cube file in current folder"
		write(*,*) "3 Export data to output.txt in current folder"
		write(*,"(a,f10.5)") " 4 Set the value of isosurface to be shown, current:",sur_value
		write(*,*) "5 Multiply all grid data by a factor"
		write(*,*) "6 Divide all grid data by a factor"
		write(*,*) "7 Add a value to all grid data"
		write(*,*) "8 Substract a value from all grid data"
		write(*,"(a)") " 9 Multiply all grid data by Hirshfeld weights of a fragment (can be used to only make isosurface around interested fragment visible)"
        if (allocated(cubmat_bk)) write(*,*) "10 Restore original grid data"
		read(*,*) i
		
		if (i==-1) then
			call drawisosurgui(1)
            
		else if (i==0) then
			exit
            
		else if (i==1) then
			idrawisosur=1
			isavepic=1
			call drawmol
			isavepic=0
			write(*,*) "Graphical file has been saved to current folder with ""dislin"" prefix"
            
		else if (i==2) then
			write(*,*) "Exporting cube file, please wait..."
            if (iaddprefix==1) call addprefix(outcubfile)
			open(10,file=outcubfile,status="replace")
			call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
			close(10)
			write(*,"(' Done! Grid data has been exported to ',a,' in current folder')") trim(outcubfile)
            if (iaddprefix==0) write(*,"(a)") " Hint: If you want to add input file name as prefix of the outputted &
            &cube file, you can set ""iaddprefix"" in settings.ini to 1"
            
		else if (i==3) then
            c200tmp="output.txt"
            if (iaddprefix==1) call addprefix(c200tmp)
			open(10,file=c200tmp,status="replace")
			write(*,"(a)") " Outputting "//trim(c200tmp)//" in current folder..."
			do k=1,nz
				do j=1,ny
					do i=1,nx
						call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
						write(10,"(3f12.6,2x,1PE18.8E3)") tmpx*b2a,tmpy*b2a,tmpz*b2a,cubmat(i,j,k)
					end do
				end do
                call showprog(k,nz)
			end do
			close(10)
			write(*,"(a)") " Output finished, column 1/2/3/4 correspond to X/Y/Z/value. The coordinate unit is Angstrom, the unit of the calculated function is a.u."
		
        else if (i==4) then
			write(*,*) "Input the value of isosurface, e.g. 0.02"
			read(*,*) sur_value
            
		else if (i==5.or.i==6.or.i==7.or.i==8) then
            allocate(cubmat_bk(nx,ny,nz))
            cubmat_bk=cubmat
			write(*,*) "Input a value, e.g. 1.5"
			read(*,*) tmpval
            if (i==5) then
				cubmat=cubmat*tmpval
            else if (i==6) then
				cubmat=cubmat/tmpval
            else if (i==7) then
				cubmat=cubmat+tmpval
            else if (i==8) then
				cubmat=cubmat-tmpval
            end if
			write(*,*) "Done! Grid data has been updated"
            
        else if (i==9) then
            allocate(cubmat_bk(nx,ny,nz))
            cubmat_bk=cubmat
			write(*,*) "Input index of the atoms in the fragment to calculate Hirshfeld weight"
			write(*,*) "e.g. 2,3,7-10,19"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,ntmp)
            if (allocated(tmparrint)) deallocate(tmparrint)
			allocate(tmparrint(ntmp))
			call str2arr(c2000tmp,ntmp,tmparrint)            
            
            write(*,*) "Calculating Hirshfeld weight at grid points..."
            iprog=0
			!$OMP PARALLEL DO SHARED(iprog,cubmat) PRIVATE(i,j,k,tmpx,tmpy,tmpz,promol,fragdens,tmpdens,weight) schedule(dynamic) NUM_THREADS(nthreads)
			do k=1,nz
				do j=1,ny
					do i=1,nx
						call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
						promol=0
                        fragdens=0
						do iatm=1,ncenter
							tmpdens=calcatmdens(iatm,tmpx,tmpy,tmpz,0)
							promol=promol+tmpdens
							if (any(tmparrint(1:ntmp)==iatm)) fragdens=fragdens+tmpdens
						end do
						if (promol==0) then
							weight=0
						else
							weight=fragdens/promol
						end if
                        cubmat(i,j,k)=cubmat(i,j,k)*weight
					end do
				end do
				!$OMP CRITICAL
                iprog=iprog+1
				call showprog(iprog,nz)
                !$OMP END CRITICAL
			end do
			!$OMP END PARALLEL DO
			write(*,*) "Done! Grid data has been updated"
            
        else if (i==10) then
			cubmat=cubmat_bk
			deallocate(cubmat_bk)
			write(*,*) "Original grid data has been restored"
		end if
	end do
end if
end subroutine





!!!----------------- Set contour lines
subroutine setcontour
use defvar
use util
implicit real*8 (a-h,o-z)
character outfilename*200,selectyn,c2000tmp*2000,c80tmp*80
integer,allocatable :: tmparr(:)
do while(.true.)
	write(*,*)
	write(*,*) "Current contour lines:"
	if (ncontour==0) then
		write(*,*) "None"
	else
		do i=1,ncontour  !The lastest contour line serial
			write(*,"(i3,':',f19.8,'   ')",advance='no') i,ctrval(i)
			if (mod(i,3)==0) write(*,*)
		end do
		if (mod(ncontour,3)/=0) write(*,*)
	end if
	if (allocated(boldlinelist)) then
		write(*,*) "Indices of bolded lines:"
		write(*,"(15i5)") boldlinelist
	end if
	write(*,*) "1 Save setting and return"
	write(*,*) "2 Change value of a contour line"
	write(*,*) "3 Add a new contour line"
	write(*,*) "4 Delete some contour lines"
    write(*,*) "5 Use built-in contour values suitable for special purpose"
	write(*,*) "6 Save contour setting to external file"
	write(*,*) "7 Load contour setting from external file"
	write(*,*) "8 Generate contour value by arithmetic progression"
	write(*,*) "9 Generate contour value by geometric series"
	if (.not.allocated(boldlinelist)) write(*,*) "10 Enable bolded style for some contour lines "
	if (allocated(boldlinelist)) write(*,*) "10 Disable bolded style"
	write(*,*) "11 Set color for positive contour lines, current: ",trim(colorname(iclrindctrpos))
	write(*,*) "12 Set line style and width for positive contour lines"
	write(*,*) "13 Set color for negative contour lines, current: ",trim(colorname(iclrindctrneg))
	write(*,*) "14 Set line style and width for negative contour lines"
	write(*,*) "15 Set line style and width suitable for publication"
	read(*,*) isel 
	if (isel==1) then
		return
	else if (isel==2) then
		write(*,*) "Input the contour line index and new value"
		write(*,*) "e.g. 4,0.015 means setting the value of the 4th contour line to 0.015"
		read(*,*) isel,selfctrval
		if (isel>=1.and.isel<=ncontour) then
			ctrval(isel)=selfctrval
		else
			write(*,*) "The number exceed valid range!"
		end if
	else if (isel==3) then
		if (ncontour<size(ctrval)) then
			write(*,*) "Input value of the contour line you want to add, e.g. 0.02"
			read(*,*) selfctrval
			ncontour=ncontour+1
			ctrval(ncontour)=selfctrval
		else
			write(*,"(a,i6)") " Error: The number of contour line could not exceed",size(ctrval)
		end if
	else if (isel==4) then
		write(*,*) "Input index of the contour lines you want to delete, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, then all contour lines will be deleted"
		read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
			ncontour=0
        else
			call str2arr(c2000tmp,ntmp)
			if (allocated(tmparr)) deallocate(tmparr)
			allocate(tmparr(ntmp))
			call str2arr(c2000tmp,ntmp,tmparr)
			if (any(tmparr<1).or.any(tmparr>ncontour)) then
				write(*,*) "Some indices you inputted exceeded valid range!"
			else
				do idx=1,ntmp
					ictr=tmparr(idx)
					ncontour=ncontour-1
					ctrval(ictr:ncontour)=ctrval(ictr+1:ncontour+1)
					where(tmparr>ictr) tmparr=tmparr-1
				end do
			end if
        end if
	else if (isel==5) then
		write(*,*) "Use which built-in contour line values?"
        write(*,*) "-1 Default contour line values"
        write(*,*) "0 Return"
        write(*,*) "1 For plotting ELF and LOL (0.0 to 1.0 with step size of 0.1)"
        write(*,*) "2 For plotting ELF and LOL (0.0 to 1.0 with step size of 0.05)"
        write(*,*) "3 For plotting orbital wavefunction (+ and - 0.01*2^(i-1), i=1-28)"
        write(*,*) "4 For plotting density difference (+ and - 0.001*2^i, i=1-20)"
        write(*,*) "5 For plotting density difference (+ and - 0.0001*2^i, i=1-23)"
        read(*,*) isel2
        if (isel2==-1) then
			call gencontour(0,0D0,0D0,0)
        else if (isel2==0) then
			cycle
        else if (isel2==1) then
			ncontour=11
            do idx=1,ncontour
				ctrval(idx)=0.1D0*(idx-1)
            end do
        else if (isel2==2) then
			ncontour=21
            do idx=1,ncontour
				ctrval(idx)=0.05D0*(idx-1)
            end do
        else if (isel2==3) then
            do idx=1,28
				ctrval(idx)=0.01D0*2**(idx-1)
				ctrval(28+idx)=-0.01D0*2**(idx-1)
            end do
			ncontour=56
        else if (isel2==4) then
            do idx=1,20
				ctrval(idx)=0.001D0*2**idx
				ctrval(20+idx)=-0.001D0*2**idx
            end do
			ncontour=40
        else if (isel2==5) then
            do idx=1,23
				ctrval(idx)=0.0001D0*2**idx
				ctrval(23+idx)=-0.0001D0*2**idx
            end do
			ncontour=46
        end if
	else if (isel==6) then
		write(*,*) "Input file path for exporting current setting, e.g. C:\ltwd.txt"
		read(*,*) outfilename
		open(10,file=outfilename)
		write(10,"(f19.8)") (ctrval(i),i=1,ncontour)
		close(10)
		write(*,"(a,a)") " Contour setting has been saved to ",trim(outfilename)
	else if (isel==7) then
		write(*,*) "Input file path, e.g. C:\ctr.txt"
		do while(.true.)
			read(*,"(a)") extctrsetting
			inquire(file=extctrsetting,exist=alive)
			if (alive) exit
			write(*,*) "File not found, input again"
		end do
		open(10,file=extctrsetting,status="old")
		ierror=0
		ncontour=0
		do i=1,size(ctrval)
			read(10,*,iostat=ierror) temp
			if (ierror/=0) exit
			ctrval(i)=temp
			ncontour=i
		end do
		close(10)
	else if (isel==8.or.isel==9) then
		write(*,*) "Input starting value, step and total number"
		if (isel==8) write(*,*) "e.g. 0.4,0.1,10 generates 0.4,0.5,0.6,0.7 ... 1.3"
		if (isel==9) write(*,*) "e.g. 2,3,10 generates 2,6,18,54 ... 39366"
		read(*,*) ctrgenstrval,ctrgenstep,igennum
		write(*,*) "If removing existing contour lines? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			ncontour=0
		else if (selectyn=='n'.or.selectyn=='N') then !Append to existed contour lines
			if (ncontour+igennum>size(ctrval)) then
				igennum=size(ctrval)-ncontour
				write(*,*) "Warning: The total number of coutour lines exceeded the upper limit!"
			end if
		end if
		do i=1,igennum
			if (isel==8) ctrval(ncontour+i)=ctrgenstrval+ctrgenstep*(i-1)
			if (isel==9) ctrval(ncontour+i)=ctrgenstrval*ctrgenstep**(i-1)
		end do
		ncontour=ncontour+igennum
	else if (isel==10) then
		if (allocated(boldlinelist)) then !Disable bold line
			deallocate(boldlinelist)
			write(*,*) "No line is bolded now, you can select this function again to set bolded lines"
		else
			write(*,*) "Set how many lines as bolded style? e.g. 5"
			read(*,*) numboldline
			allocate(boldlinelist(numboldline))
			do itmp=1,numboldline
				write(*,"(' Input the index of bolded line',i4)") itmp
				read(*,*) boldlinelist(itmp)
			end do
		end if
	else if (isel==11) then
		write(*,*) "Use which color?"
		call selcolor(iclrindctrpos)
	else if (isel==12) then
		write(*,*) "Input length of line segment and interstice"
		write(*,*) "e.g. 1,0 means solid line; 1,10 means DOT; 10,10 means DASH"
		write(*,*) "     10,15 means DASH with larger interstice"
		write(*,*) "Note: 1,0 and 10,15 are default for positive and negative lines, respectively"
        write(*,"(' Current values are:',2i5)") ctrposstyle(1),ctrposstyle(2)
        write(*,*) "If press ENTER button directly, current setting will be kept unchanged"
        read(*,"(a)") c80tmp
		if (c80tmp/=" ") read(c80tmp,*) ctrposstyle(1),ctrposstyle(2)
		write(*,*) "Input line width, e.g. 2"
        write(*,"(' Current value is',i5)") iwidthposctr
		read(*,*) iwidthposctr
	else if (isel==13) then
		write(*,*) "Use which color?"
		call selcolor(iclrindctrneg)
	else if (isel==14) then
		write(*,*) "Input length of line segment and interstice"
		write(*,*) "e.g. 1,0 means solid line; 1,10 means DOT; 10,10 means DASH"
		write(*,*) "     10,15 means DASH with larger interstice"
		write(*,*) "Note: 1,0 and 10,15 are default for positive and negative lines, respectively"
        write(*,"(' Current values are:',2i5)") ctrnegstyle(1),ctrnegstyle(2)
        write(*,*) "If press ENTER button directly, current setting will be kept unchanged"
        read(*,"(a)") c80tmp
        if (c80tmp/=" ") read(c80tmp,*) ctrnegstyle(1),ctrnegstyle(2)
		write(*,*) "Input line width, e.g. 2"
        write(*,"(' Current value is',i5)") iwidthnegctr
		read(*,*) iwidthnegctr
	else if (isel==15) then
		iclrindctrpos=11
		iclrindctrneg=3
		ctrposstyle(1)=1
		ctrposstyle(2)=0
		iwidthposctr=4
		ctrnegstyle(1)=10
		ctrnegstyle(2)=15
		iwidthnegctr=4
		write(*,"(a)") " OK! Saved picture with current line setting should be suitable for publication purpose"
	end if
end do
end subroutine



!!------------------ Set marker of critical points (CPs) and paths on contour/gradient map
subroutine settopomark
use defvar
implicit real*8 (a-h,o-z)
do while(.true.)
	write(*,*)
	write(*,*) "0 Return"
	if (imark3n3==0) write(*,*) "1 Enable showing (3,-3) CPs"
	if (imark3n3==1) write(*,*) "1 Disable showing (3,-3) CPs"
	if (imark3n1==0) write(*,*) "2 Enable showing (3,-1) CPs"
	if (imark3n1==1) write(*,*) "2 Disable showing (3,-1) CPs"
	if (imark3p1==0) write(*,*) "3 Enable showing (3,+1) CPs"
	if (imark3p1==1) write(*,*) "3 Disable showing (3,+1) CPs"
	if (imark3p3==0) write(*,*) "4 Enable showing (3,+3) CPs"
	if (imark3p3==1) write(*,*) "4 Disable showing (3,+3) CPs"
	if (imarkpath==0) write(*,*) "5 Enable showing paths"
	if (imarkpath==1) write(*,*) "5 Disable showing paths"
	write(*,"(a,f8.3,' Bohr')") " 8 Set distance threshold for showing CPs, current: ",disshowCP
	write(*,"(a,f8.3,' Bohr')") " 9 Set distance threshold for showing paths, current: ",disshowpath
	write(*,"(a,i3)") " 10 Set size of markers of CPs, current: ",sizemarkcp
	write(*,"(a,i3)") " 11 Set thickness of topology paths, current: ",sizemarkpath
	write(*,"(a,a)") " 12 Set color of topology paths, current: ",trim(colorname(iclrpath))
	write(*,"(a,i3)") " 13 Set thickness of the interbasin paths derived from (3,-1), current: ",sizemark3n1path
	write(*,"(a,a)") " 14 Set color of the interbasin paths derived from (3,-1), current: ",trim(colorname(iclr3n1path))
    write(*,*) "15 Set color for CPs"
	read(*,*) isel

	if (isel==0) then
		return
	else if (isel==1) then
		if (imark3n3==1) then
			imark3n3=0
		else
			imark3n3=1
		end if
	else if (isel==2) then
		if (imark3n1==1) then
			imark3n1=0
		else
			imark3n1=1
		end if
	else if (isel==3) then
		if (imark3p1==1) then
			imark3p1=0
		else
			imark3p1=1
		end if
	else if (isel==4) then
		if (imark3p3==1) then
			imark3p3=0
		else
			imark3p3=1
		end if
	else if (isel==5) then
		if (imarkpath==1) then
			imarkpath=0
		else
			imarkpath=1
		end if
	else if (isel==8) then
		write(*,*) "Input the distance threshold in Bohr, e.g. 0.5"
		read(*,*) disshowCP
	else if (isel==9) then
		write(*,*) "Input the distance threshold in Bohr, e.g. 0.5"
		read(*,*) disshowpath
	else if (isel==10) then
		write(*,*) "Input the size value, e.g. 30"
		read(*,*) sizemarkcp
	else if (isel==11) then
		write(*,*) "Input the thickness value, e.g. 5"
		read(*,*) sizemarkpath
	else if (isel==12) then
		write(*,*) "Use which color for topology paths?"
		call selcolor(iclrpath)
	else if (isel==13) then
		write(*,*) "Input the thickness value, e.g. 5"
		read(*,*) sizemark3n1path
	else if (isel==14) then
		write(*,*) "Use which color for topology paths?"
		call selcolor(iclr3n1path)
	else if (isel==15) then
		write(*,*) "Hint: The default CP colors can be changed by ""CP_RGB_2D"" in settings.ini"
        write(*,*)
		write(*,*) "Set color for which kind of CP?"
        write(*,*) "1 (3,-3)"
        write(*,*) "2 (3,-1)"
        write(*,*) "3 (3,+1)"
        write(*,*) "4 (3,+3)"
        read(*,*) iCPtype
        write(*,*) "Change to which color?"
        call selcolor(iclrind)
		call clridx2RGB(iclrind,Rcomp,Gcomp,Bcomp)
        if (iCPtype==1) then
			CP3n3RGB_2D=(/Rcomp,Gcomp,Bcomp/)
        else if (iCPtype==2) then
			CP3n1RGB_2D=(/Rcomp,Gcomp,Bcomp/)
        else if (iCPtype==3) then
			CP3p1RGB_2D=(/Rcomp,Gcomp,Bcomp/)
        else if (iCPtype==4) then
			CP3p3RGB_2D=(/Rcomp,Gcomp,Bcomp/)
        end if
        write(*,*) "Done!"
	end if
end do
end subroutine




!!!--------------------------  Set content of custom plot
subroutine customplotsetup
use defvar 
implicit real*8 (a-h,o-z)
write(*,*) "How many files to deal with? e.g. 3 (Excluding the already loaded file)"
read(*,*) ncustommap
if (allocated(custommapname)) deallocate(custommapname)
if (allocated(customop)) deallocate(customop)
allocate(custommapname(ncustommap))
allocate(customop(ncustommap))
write(*,*) "Avaliable operators: +,-,*,/"
write(*,"(a)") " e.g. -,C:\sob.wfn means subtracting property of C:\sob.wfn from the file loaded when Multiwfn boots up"
do i=1,ncustommap
	write(*,"(' Input operator and file path of system',i5)") i
	do while(.true.)
		read(*,"(a1,1x,a)") customop(i),custommapname(i)
		if (customop(i)/='+'.and.customop(i)/='-'.and.customop(i)/='*'.and.customop(i)/='/') then
			write(*,*) "Error: You did not properly specify an operator, input again!"
			cycle
		end if
		ltmp=len_trim(custommapname(i))
		if (custommapname(i)(1:1)=='"'.or.custommapname(i)(1:1)=="'") custommapname(i)(1:1)=" "
		if (custommapname(i)(ltmp:ltmp)=='"'.or.custommapname(i)(ltmp:ltmp)=="'") custommapname(i)(ltmp:ltmp)=" "
		inquire(file=custommapname(i),exist=alive)
		if (alive.eqv..true.) then
			exit
		else
			write(*,*) "File not found, input again"
		end if
	end do
end do
end subroutine



!!-------- Generate contour line values
!itype=0: For general case
!itype=1: Specific for ELF/LOL
!itype=2: Specify vallow and valhigh, generate evenly distributed nctr contour lines (no line corresponds to vallow, but a line corresponds to valhigh)
subroutine gencontour(itype,vallow,valhigh,nctr)
use defvar
implicit real*8 (a-h,o-z)
integer itype,nctr
real*8 vallow,valhigh
if (itype==0) then
	ncontour=62
	ctrval(1)=1D-3
	do i=0,9
		ctrval(3*i+2)=2*1D-3*10**i
		ctrval(3*i+3)=4*1D-3*10**i
		ctrval(3*i+4)=8*1D-3*10**i
	end do
	ctrval(32:62)=-ctrval(1:31)
else if (itype==1) then
	ncontour=21
	do i=1,ncontour
		ctrval(i)=(i-1)*0.05D0
	end do
else if (itype==2) then
	ncontour=nctr
    stp=(valhigh-vallow)/ncontour
	do i=1,ncontour
		ctrval(i)=vallow+i*stp
	end do
end if
end subroutine





!!----------- A general interface for plotting color-filled or contour map based on "planemat" data
!In fact this is a simplified version of post-processing menu of main function 4
!Gradient map is not available, topology information cannot be shown
!Input data:
!  mapname: The name of the map, shown as title
!  outfilename: Name of exporting file
!  xlow,xhigh,ylow,yhigh: Lower and upper limits of X and Y axis in the map
!  zlow and zhigh: Color scale range
!Before using this, below data should be initialized:
!  orgx,orgy,dx,dy: They should be defined like usual plane map
!  ngridnum1 and ngridnum2: Number of grids of planemat
!  plesel: Type of plane (XY, YZ...), like usual plane map in main function 4
!  Contour lines should be initialized by subroutine "gencontour"
!  planestpx,planestpy,planestpz: Label intervals in X,Y,Z
!  orgx2D,orgy2D,orgz2D: X, Y, Z of origin of the plane in molecular Cartesian space, used to determine if showing atomic labels, etc.
!  idrawtype and idrawcontour should be properly set, see code of this routine
subroutine planemap_interface(mapname,outfilename,xlow,xhigh,ylow,yhigh,zlow,zhigh)
use defvar
use plot
implicit real*8 (a-h,o-z)
character(len=*) mapname,outfilename
character selectyn,c80tmp*80
real*8 xlow,xhigh,ylow,yhigh,zlow,zhigh

do while(.true.)
    write(*,*)
    ntmp=(63-len(mapname))/2
    write(*,"(' ')",advance='no')
    do itmp=1,ntmp
        write(*,"('-')",advance='no')
    end do
    write(*,"(' Plotting ',a,' map ')",advance='no') mapname
    do itmp=1,ntmp
        write(*,"('-')",advance='no')
    end do
    write(*,*) "-5 Save (load) all plotting settings to (from) from an external file"
    write(*,*) "-4 Multiply the current data by a factor"
    write(*,*) "-3 Change other plotting settings"
    write(*,*) "-2 Export plane data as "//outfilename//".txt in current folder"
    write(*,*) "-1 Return"
    write(*,*) "0 Show the map on screen"
    write(*,*) "1 Save the map as graphical file in current folder"
    if (idrawtype==1.and.idrawcontour==0) then
        write(*,*) "2 Choose map type, current: Color-filled map"
    else if (idrawtype==1.and.idrawcontour==1) then
        write(*,*) "2 Choose map type, current: Color-filled map with contour lines"
    else if (idrawtype==2) then
        write(*,*) "2 Choose map type, current: Contour line map"
    end if
	if (iatom_on_plane==0) write(*,*) "3 Toggle showing atom labels, current: No"
	if (iatom_on_plane==1) write(*,*) "3 Toggle showing atom labels, current: Yes"
    if (iatom_on_plane==1) write(*,"(a,f8.3,' Bohr')") " 4 Set distance threshold for showing atom labels, current:",disshowlabel
    if (ibond_on_plane==0) write(*,*) "5 Toggle showing bonds, current: No"
	if (ibond_on_plane==1) write(*,*) "5 Toggle showing bonds, current: Yes"
	if (ilenunit2D==1) write(*,*) "6 Change length unit of the graph to Angstrom"
	if (ilenunit2D==2) write(*,*) "6 Change length unit of the graph to Bohr"
	if (idrawtype==1) then
		write(*,"(a,2f7.3,f12.6)") " 7 Set label intervals in X, Y and color scale axes, current:",planestpx,planestpy,planestpz
	else
		write(*,"(a,2f7.3)") " 7 Set label intervals in X and Y axes, current:",planestpx,planestpy
	end if
    !Options specific for color-filled map
    if (idrawtype==1) then
		if (abs(zlow)<1000000.and.abs(zhigh)<1000000) then
			write(*,"(a,2f15.7)") " 8 Set lower&upper limit of color scale, current:",zlow,zhigh
		else
			write(*,"(a,2(1PE15.6))") " 8 Set lower&upper limit of color scale, current:",zlow,zhigh
		end if
        write(*,"(a,a)") " 9 Set color transition, current: ",trim(clrtransname(iclrtrans))
    end if
    if (idrawcontour==1) then
		write(*,*) "10 Change contour line setting"
		if (ilabel_on_contour==0) write(*,*) "11 Toggle showing isovalue on contour lines, current: No"
		if (ilabel_on_contour==1) write(*,*) "11 Toggle showing isovalue on contour lines, current: Yes"
    end if
    read(*,*) isel
    
    if (isel==-5) then
        call saveload2Dplottingsetting(zlow,zhigh)
	else if (isel==-4) then		
		write(*,*) "Input the value to be multiplied to the current plane data, e.g. 0.3"
		read(*,*) scaleval
        planemat=planemat*scaleval
    else if (isel==-3) then
        if (idrawtype==1) then
            call plane_othersetting(1) !Involve Z-axis
        else if (idrawtype==2) then
            call plane_othersetting(0) !Do not involve Z-axis
        end if
    else if (isel==-2) then
        open(10,file=outfilename//".txt",status="replace")
        do ix=1,nx
            xpos=orgx2D+(ix-1)*dx
            do iy=1,ny
                ypos=orgy2D+(iy-1)*dy
                write(10,"(2f12.6,f16.8)") xpos*b2a,ypos*b2a,planemat(ix,iy)
            end do
        end do
        close(10)
        write(*,"(a)") " Done! plane data has been exported to "//outfilename//".txt in current folder. &
        &The first two column correspond to X and Y coordinates in Angstrom, the last column is data"
    else if (isel==-1) then
        return
    else if (isel==0.or.isel==1) then
        if (isel==0) then
            isavepic=0
        else if (isel==1) then
            isavepic=1
        end if
		call drawplane(xlow,xhigh,ylow,yhigh,zlow,zhigh)
        if (isel==1) then
            isavepic=0
            write(*,"(a,a,a)") " The map have been saved as ",trim(graphformat)," format with ""dislin"" prefix in current folder"
        end if
    else if (isel==2) then
        write(*,*) "Choose map type"
        write(*,*) "1 Color-filled map"
        write(*,*) "2 Color-filled map with contour lines"
        write(*,*) "3 Contour line map"
        read(*,*) isel2
        if (isel2==1.or.isel2==2) then
            idrawtype=1
            if (isel2==1) then
                idrawcontour=0
            else if (isel2==2) then
                idrawcontour=1
            end if
        else if (isel2==3) then
            idrawtype=2
            idrawcontour=1
        end if
    else if (isel==3) then
		if (iatom_on_plane==1) then
			iatom_on_plane=0
		else if (iatom_on_plane==0) then
			iatom_on_plane=1
			write(*,*) "Use which color for labelling atoms?"
			call selcolor(iclrindatmlab)
		end if
    else if (isel==4) then
		write(*,*) "Input distance threshold for plotting atomic labels (in Bohr), e.g. 0.5"
        write(*,*) "If you want to input in Angstrom, add ""A"" suffix, e.g. 0.45 A"
		write(*,*) "Note: The default value can be set by ""disshowlabel"" in settings.ini"
        read(*,"(a)") c80tmp
		read(c80tmp,*) disshowlabel
        if (index(c80tmp,'A')/=0.or.index(c80tmp,'a')/=0) disshowlabel=disshowlabel/b2a
		write(*,"(a)") " If also show labels of the atoms that beyond this criterion as light face type? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			iatom_on_plane_far=1
		else
			iatom_on_plane_far=0
		end if
    else if (isel==5) then
		if (ibond_on_plane==0) then
			ibond_on_plane=1
			write(*,"(a,/)") " Note: The bonding will be empirically determined according to interatomic distance and atomic covalent radii"
			write(*,*) "Use which color for drawing the bonds?"
			call selcolor(iclrindbndlab)
		else
			ibond_on_plane=0
		end if
    else if (isel==6) then
		if (ilenunit2D==1) then
			ilenunit2D=2
		else if (ilenunit2D==2) then
			ilenunit2D=1
		end if
    else if (isel==7) then
		if (idrawtype==1) then
			write(*,"(a)") " Input interval between the labels in X, Y and color scale axes, e.g. 1.5,2.0,0.1"
			read(*,*) planestpx,planestpy,planestpz
		else if (idrawtype==2) then
			write(*,"(a)") " Input interval between the labels in X and Y axes, e.g. 1.5,2.0"
			read(*,*) planestpx,planestpy
		end if
    else if (isel==8) then
		write(*,*) "Input lower & upper limit of color scale, e.g. -0.3,0.3"
		read(*,*) zlow,zhigh
    else if (isel==9) then
        call selcolortable
    else if (isel==10) then
        call setcontour
    else if (isel==11) then
		if (ilabel_on_contour==1) then
			ilabel_on_contour=0
		else if (ilabel_on_contour==0) then
			ilabel_on_contour=1
			write(*,*) "Input label size, e.g. 30"
			read(*,*) ictrlabsize
			write(*,"(a)") " Hint: The number of digits after the decimal point of label on contour lines can be set by ""numdigctr"" in settings.ini"
		end if
    end if
end do      
end subroutine



!!---------- Save/load plotting settings to/from an .txt file
!The clrlow and clrhigh are lower and upper limits of color scale. They are arguments since they are not global variables
subroutine saveload2Dplottingsetting(clrlow,clrhigh)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 clrlow,clrhigh
character c200tmp*200,strtmp*20

do while(.true.)
    write(*,*)
    write(*,*) "0 Return"
    write(*,*) "1 Save plotting settings to specific .txt file in current folder"
    write(*,*) "2 Load plotting settings from specific .txt file in current folder"
    read(*,*) isel

    if (isel==0) then
        return
    else if (isel==1) then
        write(*,"(a)") " Input the path of the file to save plotting settings, for example, C:\RAS\masking.txt"
        write(*,*) "If pressing ENTER button directly, the file will be planeplot.txt"
        read(*,"(a)") c200tmp
        if (c200tmp==" ") c200tmp="planeplot.txt"
        open(10,file=c200tmp,status="replace")
        !write(10,"('idrawtype',i5)") idrawtype
        !write(10,"('graph2Dwidth',i5)") graph2Dwidth
        !write(10,"('graph2Dheight',i5)") graph2Dheight
        !write(10,"('graphformat ',a)") graphformat
        write(10,"('clrlow',1PE18.8)") clrlow
        write(10,"('clrhigh',1PE18.8)") clrhigh
        write(10,"('iclrtrans',i5)") iclrtrans
        write(10,"('idrawcontour',i5)") idrawcontour
        write(10,"('iatom_on_plane',i5)") iatom_on_plane
        write(10,"('disshowlabel',f12.6)") disshowlabel
        write(10,"('iatom_on_plane_far',i5)") iatom_on_plane_far
        write(10,"('ibond_on_plane',i5)") ibond_on_plane
        write(10,"('ilenunit2D',i5)") ilenunit2D
        write(10,"('planestpx',f12.6)") planestpx
        write(10,"('planestpy',f12.6)") planestpy
        write(10,"('planestpz',f12.6)") planestpz
        write(10,"('ilabel_on_contour',i5)") ilabel_on_contour
        write(10,"('iclrindatmlab',i6)") iclrindatmlab
        write(10,"('iclrindbndlab',i6)") iclrindbndlab
        !Other settings
        write(10,"('itickreverse',i5)") itickreverse
        write(10,"('numdigx',i5)") numdigx
        write(10,"('numdigy',i5)") numdigy
        write(10,"('numdigz',i5)") numdigz
        write(10,"('iticks',i5)") iticks
        write(10,"('plane_axistextsize',i5)") plane_axistextsize
        write(10,"('plane_axisnamesize',i5)") plane_axisnamesize
        write(10,"('iatmlabtype',i5)") iatmlabtype
        write(10,"('pleatmlabsize',i5)") pleatmlabsize
        write(10,"('bondthick2D',i5)") bondthick2D
        !Related to relief map
        write(10,"('surcolorzmin',f16.8)") surcolorzmin
        write(10,"('surcolorzmax',f16.8)") surcolorzmax
        write(10,"('drawsurmesh ',a)") drawsurmesh
        !Related to gradient line map
        write(10,"('igrad_arrow',i6)") igrad_arrow
        write(10,"('gradplotstep',f8.4)") gradplotstep
        write(10,"('gradplotdis',f8.4)") gradplotdis
        write(10,"('gradplottest',f8.4)") gradplottest
        write(10,"('iclrindgradline',i6)") iclrindgradline
        write(10,"('iwidthgradline',i6)") iwidthgradline
        write(10,"('stream_intmethod ',a)") stream_intmethod
        !Related to vector field map
        write(10,"('cutgradvec',f8.4)") cutgradvec
        write(10,"('igrad_arrow',i6)") igrad_arrow
        write(10,"('icolorvecfield',i6)") icolorvecfield
        write(10,"('iinvgradvec',f8.4)") iinvgradvec
        write(10,"('vecclrind',f8.4)") vecclrind
        !Related to plotting topology information
        write(10,"('imark3n3',i6)") imark3n3
        write(10,"('imark3n1',i6)") imark3n1
        write(10,"('imark3p1',i6)") imark3p1
        write(10,"('imark3p3',i6)") imark3p3
        write(10,"('imarkpath',i6)") imarkpath
        write(10,"('disshowCP',f10.3)") disshowCP
        write(10,"('disshowpath',f10.3)") disshowpath
        write(10,"('sizemarkcp',i6)") sizemarkcp
        write(10,"('sizemarkpath',i6)") sizemarkpath
        write(10,"('iclrpath',i6)") iclrpath
        write(10,"('sizemark3n1path',i6)") sizemark3n1path
        write(10,"('iclr3n1path',i6)") iclr3n1path
        write(10,"('CP3n3RGB_2D',3f10.6)") CP3n3RGB_2D
        write(10,"('CP3n1RGB_2D',3f10.6)") CP3n1RGB_2D
        write(10,"('CP3p1RGB_2D',3f10.6)") CP3p1RGB_2D
        write(10,"('CP3p3RGB_2D',3f10.6)") CP3p3RGB_2D
        !Related to contour lines
        write(10,"('ictrlabsize',i5)") ictrlabsize
        write(10,"('iclrindctrpos',i5)") iclrindctrpos
        write(10,"('iclrindctrneg',i5)") iclrindctrneg
        write(10,"('iwidthposctr',i5)") iwidthposctr
        write(10,"('iwidthnegctr',i5)") iwidthnegctr
        write(10,"('ctrposstyle',2i5)") ctrposstyle(:)
        write(10,"('ctrnegstyle',2i5)") ctrnegstyle(:)
        write(10,"('ifillctrline',i5)") ifillctrline
        write(10,"('ishowclrfill_bar',i5)") ishowclrfill_bar
        write(10,"('ncontour',i5)") ncontour
        do ictr=1,ncontour
	        write(10,"(1PE18.8)") ctrval(ictr)
        end do
        close(10)
        write(*,"(a)") " Done! Plotting settings have been saved to "//trim(c200tmp)
    else if (isel==2) then
        write(*,"(a)") " Input the path of the file from which the plotting settings will be loaded, e.g. C:\RAS\masking.txt"
        write(*,*) "If pressing ENTER button directly, the file will be planeplot.txt"
        do while(.true.)
            read(*,"(a)") c200tmp
            if (c200tmp==" ") c200tmp="planeplot.txt"
            inquire(file=c200tmp,exist=alive)
	        if (alive) exit
	        write(*,*) "Cannot find the file, input again!"
        end do
        open(10,file=c200tmp,status="old")
        read(10,*) strtmp,clrlow
        read(10,*) strtmp,clrhigh
        read(10,*) strtmp,iclrtrans
        read(10,*) strtmp,idrawcontour
        read(10,*) strtmp,iatom_on_plane
        read(10,*) strtmp,disshowlabel
        read(10,*) strtmp,iatom_on_plane_far
        read(10,*) strtmp,ibond_on_plane
        read(10,*) strtmp,ilenunit2D
        read(10,*) strtmp,planestpx
        read(10,*) strtmp,planestpy
        read(10,*) strtmp,planestpz
        read(10,*) strtmp,ilabel_on_contour
        read(10,*) strtmp,iclrindatmlab
        read(10,*) strtmp,iclrindbndlab
        !Other settings
        read(10,*) strtmp,itickreverse
        read(10,*) strtmp,numdigx
        read(10,*) strtmp,numdigy
        read(10,*) strtmp,numdigz
        read(10,*) strtmp,iticks
        read(10,*) strtmp,plane_axistextsize
        read(10,*) strtmp,plane_axisnamesize
        read(10,*) strtmp,iatmlabtype
        read(10,*) strtmp,pleatmlabsize
        read(10,*) strtmp,bondthick2D
        !Related to relief map
        read(10,*) strtmp,surcolorzmin
        read(10,*) strtmp,surcolorzmax
        read(10,*) strtmp,drawsurmesh
        !Related to gradient line map
        read(10,*) strtmp,igrad_arrow
        read(10,*) strtmp,gradplotstep
        read(10,*) strtmp,gradplotdis
        read(10,*) strtmp,gradplottest
        read(10,*) strtmp,iclrindgradline
        read(10,*) strtmp,iwidthgradline
        read(10,*) strtmp,stream_intmethod
        !Related to vector field map
        read(10,*) strtmp,cutgradvec
        read(10,*) strtmp,igrad_arrow
        read(10,*) strtmp,icolorvecfield
        read(10,*) strtmp,iinvgradvec
        read(10,*) strtmp,vecclrind
        !Related to plotting topology information
        read(10,*) strtmp,imark3n3
        read(10,*) strtmp,imark3n1
        read(10,*) strtmp,imark3p1
        read(10,*) strtmp,imark3p3
        read(10,*) strtmp,imarkpath
        read(10,*) strtmp,disshowCP
        read(10,*) strtmp,disshowpath
        read(10,*) strtmp,sizemarkcp
        read(10,*) strtmp,sizemarkpath
        read(10,*) strtmp,iclrpath
        read(10,*) strtmp,sizemark3n1path
        read(10,*) strtmp,iclr3n1path
        call readoption_vec_float(10,"CP3n3RGB_2D",' ',CP3n3RGB_2D)
        call readoption_vec_float(10,"CP3n1RGB_2D",' ',CP3n1RGB_2D)
        call readoption_vec_float(10,"CP3p1RGB_2D",' ',CP3p1RGB_2D)
        call readoption_vec_float(10,"CP3p3RGB_2D",' ',CP3p3RGB_2D)
        !Related to contour lines
        call readoption_int(10,"ictrlabsize",' ',ictrlabsize)
        call readoption_int(10,"iclrindctrpos",' ',iclrindctrpos)
        call readoption_int(10,"iclrindctrneg",' ',iclrindctrneg)
        call readoption_int(10,"iwidthposctr",' ',iwidthposctr)
        call readoption_int(10,"iwidthnegctr",' ',iwidthnegctr)
        call readoption_vec_int(10,"ctrposstyle",' ',ctrposstyle)
        call readoption_vec_int(10,"ctrnegstyle",' ',ctrnegstyle)
        call readoption_int(10,"ifillctrline",' ',ifillctrline)
        call readoption_int(10,"ishowclrfill_bar",' ',ishowclrfill_bar)
        call readoption_int(10,"ncontour",' ',ncontour)
	    read(10,*) ctrval(1:ncontour)
        close(10)
        write(*,"(a)") " Done! Plotting settings have been retrieved from "//trim(c200tmp)
    end if
end do
end subroutine