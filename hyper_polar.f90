!!!-------------- Various (hyper)polarizability analyses   
subroutine hyper_polarizability
do while(.true.)
    write(*,*)
    write(*,*) "     -------------- Various (hyper)polarizability analyses --------------"
    write(*,*) "0 Return to main menu"
	write(*,"(a)") " 1 Parse output file of (hyper)polarizability task of Gaussian and calculate various related quantities"
    write(*,"(a)") " 2 Study (hyper)polarizability by sum-over-states (SOS) method and perform two/three-level analysis"
    write(*,*) "3 (hyper)polarizability density analysis"
    write(*,*) "4 Calculate atomic polarizabilities in molecule"
    write(*,*) "5 Visualize (hyper)polarizability via unit sphere and vector representations"
    read(*,*) isel
    if (isel==0) then
        exit
    else if (isel==1) then
        call parseGauPolar
    else if (isel==2) then
        call SOS
    else if (isel==3) then
        call hyper_polar_dens
    else if (isel==4) then
        write(*,"(a)") " Please use subfunction 13 of fuzzy analysis module (main function 15) to conduct this analysis. &
        See Section 3.18.12 of manual for introduction and Section 4.15.4 for example"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
    else if (isel==5) then
        call vis_hypol
    end if
end do
end subroutine










!!-----------------------------------------------------------------------------
!!--------- Parse the output of (hyper)polarizability task of Gaussian --------
!!-----------------------------------------------------------------------------
subroutine parseGauPolar
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 dipole(3),alpha(3,3),beta(3,3,3),gamma(3,3,3,3)
real*8 alpha_org(3,3),beta_org(3,3,3),gamma_org(3,3,3,3) !The (hyper)polarizability tensor may be modified during processing, store the original ones in a.u.
real*8 eigvecmat(3,3),eigval(3),freqval(100000)
character :: c200tmp*200,sepchar,c210tmp*210,selectyn,lb(3)=(/ "X","Y","Z" /)
character(len=20) :: form,formau="(a,f17.5)",formother="(a,1PE16.6)"
integer :: irdfreq=0,ides=6,iunit=1,ioutpol=0
write(*,*) "Note: This function only works for ""polar"" tasks of Gaussian with #P"
10 continue
alpha=0D0
beta=0D0
gamma=0D0
do while(.true.)
    write(*,*)
    write(*,*) "------------ Parsing (hyper)polarizability calculated by Gaussian ------------"
	write(*,"(a,/)") " Select actual case of your Gaussian task. Option 2,4,6 print alpha, &
    option 1,3,5 print both alpha and beta, option 7 prints both alpha and gamma. -1 can be chosen &
    only if CPHF=RdFreq or polar=DCSHG was used in calculation"
    if (ioutpol==1) write(*,*) "-4 Export (hyper)polarizability as .txt file after parsing, current: Yes"
    if (ioutpol==0) write(*,*) "-4 Export (hyper)polarizability as .txt file after parsing, current: No"
	if (iunit==1) write(*,*) "-3 Set unit in the output, current: a.u."
	if (iunit==2) write(*,*) "-3 Set unit in the output, current: SI"
	if (iunit==3) write(*,*) "-3 Set unit in the output, current: esu"
	if (ides==6) write(*,*) "-2 Set output destination, current: Output to screen"
	if (ides==11) write(*,*) "-2 Set output destination, current: polar.txt in current folder"
	if (irdfreq==1) write(*,*) "-1 Toggle loading frequency-dependent result for options 1 and 7, current: Yes"
	if (irdfreq==0) write(*,*) "-1 Toggle loading frequency-dependent result for options 1 and 7, current: No"
	write(*,*) "0 Return"
	write(*,*) "1 ""Polar"" + analytic 3-order deriv. (HF,pure/hybrid DFT,semi-empirical)"
	write(*,*) "2 ""Polar"" + analytic 2-order deriv. (MP2,double-hybrid DFT,TDDFT,CIS...)"
	write(*,*) "3 ""Polar=Cubic"" + analytic 2-order deriv."
	write(*,*) "4 ""Polar"" + analytic 1-order deriv. (CISD,QCISD,CCSD,MP3,MP4(SDQ)...)"
	write(*,*) "5 ""Polar=DoubleNumer"" or ""Polar=EnOnly"" + analytic 1-order deriv."
	write(*,*) "6 ""Polar"" + energy only (CCSD(T),QCISD(T),MP4(SDTQ),MP5...)"
	write(*,*) "7 ""Polar=gamma"" + analytic 3-order deriv. (HF,DFT,semi-empirical)"
	read(*,*) isel
	if (isel==-1) then
		if (irdfreq==1) then
			irdfreq=0
		else
			irdfreq=1
		end if
	else if (isel==-2) then
		write(*,*) "6: Output to screen"
		write(*,*) "11: Output to polar.txt in current folder"
		read(*,*) ides
	else if (isel==-3) then
		write(*,*) "1: Atomic unit"
		write(*,*) "2: SI unit (C^2*m^2/J for alpha, C^3*m^3/J for beta)"
		write(*,*) "3: esu"
		read(*,*) iunit
	else if (isel==-4) then
        if (ioutpol==1) then
            ioutpol=0
        else
            ioutpol=1
        end if
	else if (isel==0) then
		return
	else
		exit
	end if
end do

if (irdfreq==1.and.(isel/=1.and.isel/=7)) then
	write(*,*) "Error: Frequency-dependent values are only available for HF/DFT/semi-empirical!"
	return
end if

open(10,file=filename,status="old")
if (ides==11) open(ides,file="polar.txt",status="replace")

if (iunit==1) then
	write(ides,*) "Note: All data shown below are in a.u."
	form=formau
else if (iunit==2) then
	write(ides,"(a)") " Note: All data shown below are in SI unit (C^2*m^2/J for alpha, C^3*m^3/J for beta)"
	form=formother
else if (iunit==3) then
	write(ides,*) "Note: All data shown below are in esu unit"
	form=formother
end if
write(ides,*)

! Dipole moment part (miu)
call loclabel(10,"Dipole moment (field-independent basis, Debye)",ifound)
if (ifound==1) then
	read(10,*)
	read(10,*) c200tmp,xtmp,c200tmp,ytmp,c200tmp,ztmp
	dipole(1)=xtmp/au2debye !Convert to a.u.
	dipole(2)=ytmp/au2debye
	dipole(3)=ztmp/au2debye
	if (iunit==2) dipole=dipole*8.47835D-30
	if (iunit==3) dipole=dipole*2.54175D-18	
	dipnorm=dsqrt(sum(dipole**2))
	write(ides,*) "Dipole moment:"
	if (iunit==1) then
		write(ides,"(' X,Y,Z=',3f12.6,'   Norm=',f12.6)") dipole(:),dipnorm
	else
		write(ides,"(' X,Y,Z=',3(1PE15.6),'   Norm=',1PE15.6)") dipole(:),dipnorm
	end if
	write(ides,*)
end if

!Selecting the result at which frequency will be loaded
if (irdfreq==1) then
	nfreqval=0
	rewind(10)
	do while(.true.)
		call loclabel(10,"-- Alpha(-w,w) frequency",ifound,0) !Not rewind
		if (ifound==1) then
			read(10,"(46x,f12.6)") tmpval
			if ( any(freqval(1:nfreqval)==tmpval) ) exit !The output in stage 2 and stage 3 are identical, so determine if has entered stage 3
			nfreqval=nfreqval+1
			freqval(nfreqval)=tmpval
		else
			exit
		end if
	end do
	write(*,*) "Load which one? Input the index, e.g. 2"
	do i=1,nfreqval
		if (freqval(i)==0) then
            write(*,"(i8,'   w=',f12.6,' (     Static    )')") i,freqval(i)
		else if (freqval(i)>0) then
            write(*,"(i8,'   w=',f12.6,' (',f12.2,'nm )')") i,freqval(i),1239.842D0/(freqval(i)*au2eV)
        end if
	end do
	read(*,*) ifreq
	if (freqval(ifreq)>0) write(*,"(' Note: Printed (hyper)polarizability will correspond to w=',f12.6,' (',f12.2,'nm )',/)") freqval(ifreq),1239.842D0/(freqval(ifreq)*au2eV)
	if (freqval(ifreq)==0) write(*,"(' Note: Printed (hyper)polarizability will correspond to w=',f12.6,' ( Static )',/)") freqval(ifreq)
	rewind(10) !Move to the beginning
end if


!!!! Print polarizability part (Alpha)
ireadalpha_archive=0
if (isel==1.or.isel==4.or.isel==5.or.isel==7) then !Standard orientation
	if (isel==1.or.isel==7) then
		if (irdfreq==0) then 
			call loclabel(10,"SCF Polarizability",ifound)
            if (ifound==0) then
                write(*,"(a)") " Error: Unable to find needed information! Please double check your keywords, and make sure that #P has been specified"
                write(*,*) "Press ENTER button to return"
                read(*,*)
                return
            end if
		else if (irdfreq==1) then
			do i=1,ifreq
				call loclabel(10,"SCF Polarizability",ifound,0) !Not rewind
                if (ifound==0) then
                    write(*,"(a)") " Error: Unable to find needed information! Please double check your keywords, and make sure that #P has been specified"
                    write(*,*) "Press ENTER button to return"
                    read(*,*)
                    return
                end if
				read(10,*)
			end do
			backspace(10)
		end if
	else if (isel==4.or.isel==5) then
		call loclabel(10,"Isotropic polarizability",ifound)
	end if
	call skiplines(10,2)
	read(10,*) rnouse,alpha(1,1)
	read(10,*) rnouse,alpha(2,1),alpha(2,2)
	read(10,*) rnouse,alpha(3,1),alpha(3,2),alpha(3,3)
else if (isel==2.or.isel==3) then !Standard orientation
	call loclabel(10,"Exact polarizability")
    read(10,"(a)") c200tmp
	read(c200tmp(24:),*,iostat=ierror) alpha(1,1),alpha(2,1),alpha(2,2),alpha(3,1),alpha(3,2),alpha(3,3)
    !Sometimes magnitude of polarizability is quite large, leading to *** in this part. In this case load data from archive part instead
    if (ierror/=0) ireadalpha_archive=1 
end if
if (isel==6.or.ireadalpha_archive==1) then !Find result from archive part, this corresponds to input orientation
    call loclabel(10,"1\1\",ifound) !First assume to be Linux output
    if (ifound==1) then
        sepchar="\"
    else
        sepchar="|"
        call loclabel(10,"1|1|")
    end if
	do while(.true.)
		read(10,"(1x,a70)") c210tmp(1:70) !Combine two lines
		read(10,"(1x,a70)") c210tmp(71:140)
		read(10,"(1x,a70)") c210tmp(141:210)
		if (index(c210tmp(1:140),sepchar//"Polar=")/=0) exit
		backspace(10)
		backspace(10)
	end do
	do istart=1,204
		if (c210tmp(istart:istart+6)==sepchar//"Polar=") exit
	end do
	do iend=istart+1,210
		if (c210tmp(iend:iend)==sepchar) exit
	end do
	c210tmp(1:istart+6)=" " !Clean other information so that the data can be read in free format
	c210tmp(iend:)=" "
	read(c210tmp,*) alpha(1,1),alpha(2,1),alpha(2,2),alpha(3,1),alpha(3,2),alpha(3,3)
end if

alpha(1,2)=alpha(2,1)
alpha(1,3)=alpha(3,1)
alpha(2,3)=alpha(3,2)
alpha_org=alpha
!Convert to other unit
if (iunit==2) alpha=alpha*1.6488D-41
if (iunit==3) alpha=alpha*1.4819D-25
if (irdfreq==0) write(ides,*) "Static polarizability:"
if (irdfreq==1) write(ides,*) "Frequency-dependent polarizability:"
write(ides,form) " XX=",alpha(1,1)
write(ides,form) " XY=",alpha(1,2)
write(ides,form) " YY=",alpha(2,2)
write(ides,form) " XZ=",alpha(1,3)
write(ides,form) " YZ=",alpha(2,3)
write(ides,form) " ZZ=",alpha(3,3)
alphaiso=(alpha(1,1)+alpha(2,2)+alpha(3,3))/3D0
write(ides,form) ' X component of polarizability:',sum(alpha(1,:))
write(ides,form) ' Y component of polarizability:',sum(alpha(2,:))
write(ides,form) ' Z component of polarizability:',sum(alpha(3,:))
write(ides,form) ' Isotropic average polarizability:',alphaiso
write(ides,"(' Isotropic average polarizability volume:',f15.6,' Angstrom^3')") alphaiso*0.14818470D0
term1=(alpha(1,1)-alpha(2,2))**2 + (alpha(1,1)-alpha(3,3))**2 + (alpha(2,2)-alpha(3,3))**2
term2=6*(alpha(1,2)**2+alpha(1,3)**2+alpha(2,3)**2)
alphaani1=dsqrt((term1+term2)/2D0)
write(ides,form) ' Polarizability anisotropy (definition 1):',alphaani1
alphaani2=dsqrt(term1/2D0)
write(ides,form) ' Polarizability anisotropy (definition 2):',alphaani2
call diagmat(alpha,eigvecmat,eigval,300,1D-10)
call sort(eigval)
if (iunit==1) then
	write(ides,"(a,3f13.4)") ' Eigenvalues of polarizability tensor:',eigval
else
	write(ides,"(a,3(1PE13.5))") ' Eigenvalues of polarizability tensor:',eigval
end if
alphaani3=eigval(3)-(eigval(1)+eigval(2))/2D0
write(ides,form) ' Polarizability anisotropy (definition 3):',alphaani3
write(ides,*)
if (isel==6.or.ireadalpha_archive==1) then
    write(ides,*) "Note: The polarizability printed above corresponds to input orientation"
else
    write(ides,*) "Note: The polarizability printed above corresponds to standard orientation"
end if
write(ides,*)
if (ioutpol==1) then
    open(12,file="alpha.txt",status="replace")
    write(12,"(1PE18.8E3)") ((alpha_org(i,j),j=1,3),i=1,3)
    close(12)
    write(*,*) "Alpha tensor has been exported to alpha.txt in current folder"
    write(*,"(a)") " The sequence is ((alpha(i,j),j=1,3),i=1,3), where j loops first. The unit is a.u."
    write(*,*)
end if

!!!! Print first hyperpolarizability part (Beta)
if (isel==1.or.isel==3.or.isel==5) then
	if (irdfreq==0) then
		if (isel==1) then
			call loclabel(10,"SCF Static Hyperpolarizability",ifound)
		else if (isel==3) then
			call loclabel(10,"Final packed hyperpolarizability",ifound)
		else if (isel==5) then
			call loclabel(10,"Static Hyperpolarizability",ifound)
		end if
		call skiplines(10,3)
		read(10,*) rnouse,beta(1,1,1) !XXX
		call skiplines(10,2)
		read(10,*) rnouse,beta(1,1,2) !XXY
		read(10,*) rnouse,beta(1,2,2),beta(2,2,2) !XYY,YYY
		call skiplines(10,2)
		read(10,*) rnouse,beta(1,1,3) !XXZ
		read(10,*) rnouse,beta(1,2,3),beta(2,2,3) !XYZ,YYZ
		read(10,*) rnouse,beta(1,3,3),beta(2,3,3),beta(3,3,3) !XZZ,YZZ,ZZZ
		!XYX=YXX    =XXY
		beta(1,2,1)=beta(1,1,2)
		beta(2,1,1)=beta(1,1,2)
		!YXY=YYX    =XYY
		beta(2,1,2)=beta(1,2,2)
		beta(2,2,1)=beta(1,2,2)
		!ZXX=XZX    =XXZ
		beta(3,1,1)=beta(1,1,3)
		beta(1,3,1)=beta(1,1,3)
		!XZY=YXZ=ZYX=YXZ=YZX   =XYZ
		beta(1,3,2)=beta(1,2,3)
		beta(2,1,3)=beta(1,2,3)
		beta(3,2,1)=beta(1,2,3)
		beta(2,1,3)=beta(1,2,3)
		beta(2,3,1)=beta(1,2,3)
		!ZYY=YZY    =YYZ
		beta(3,2,2)=beta(2,2,3)
		beta(2,3,2)=beta(2,2,3)
		!ZXZ,ZZX    =XZZ
		beta(3,1,3)=beta(1,3,3)
		beta(3,3,1)=beta(1,3,3)
		!ZYZ=ZZY    =YZZ
		beta(3,2,3)=beta(2,3,3)
		beta(3,3,2)=beta(2,3,3)
		write(ides,"(a,/)") " Note: It is well known that the sign of hyperpolarizability of Gaussian should be inverted, the outputs shown below have already been corrected"
		beta=-beta
        beta_org=beta
		!Convert to other unit
		if (iunit==2) beta=beta*3.20636D-53
		if (iunit==3) beta=beta*8.63922D-33
		write(ides,*) "Static first hyperpolarizability:"
		write(ides,form) " XXX=",beta(1,1,1)
		write(ides,form) " XXY=",beta(1,1,2)
		write(ides,form) " XYY=",beta(1,2,2)
		write(ides,form) " YYY=",beta(2,2,2)
		write(ides,form) " XXZ=",beta(1,1,3)
		write(ides,form) " XYZ=",beta(1,2,3)
		write(ides,form) " YYZ=",beta(2,2,3)
		write(ides,form) " XZZ=",beta(1,3,3)
		write(ides,form) " YZZ=",beta(2,3,3)
		write(ides,form) " ZZZ=",beta(3,3,3)
		write(ides,*)
		betaX=beta(1,1,1)+beta(1,2,2)+beta(1,3,3)
		betaY=beta(2,2,2)+beta(2,1,1)+beta(2,3,3)
		betaZ=beta(3,3,3)+beta(3,2,2)+beta(3,1,1)
		if (iunit==1) then
			write(ides,"(' Beta_X=',f16.5,'  Beta_Y=',f16.5,'  Beta_Z=',f16.5)") betaX,betaY,betaZ
		else
			write(ides,"(' Beta_X=',1PE15.5,'  Beta_Y=',1PE15.5,'  Beta_Z=',1PE15.5)") betaX,betaY,betaZ
		end if
		write(ides,form) " Magnitude of first hyperpolarizability:",dsqrt(betaX**2+betaY**2+betaZ**2)
		betaprj=(betaX*dipole(1)+betaY*dipole(2)+betaZ*dipole(3))/dipnorm
		write(ides,form) " Projection of beta on dipole moment:",betaprj
		write(ides,form) " Beta ||     :",betaprj/5D0*3D0
		write(ides,form) " Beta ||(z)  :",betaZ/5D0*3D0
		write(ides,form) " Beta _|_(z) :",betaZ/5D0
        write(ides,"(/,a)") " Note: The beta information printed above corresponds to standard orientation"
		
	else if (irdfreq==1) then !Frequency-dependent hyperpolarizability, only available for HF/DFT/semi-empirical
		write(*,*) "Loading which type of hyperpolarizability?"
		write(*,*) "1: Beta(-w;w,0)   2: Beta(-2w;w,w)"
		write(*,*) "Note: Option 2 is meaningless if ""DCSHG"" or gamma keyword was not used" 
		read(*,*) ibeta
		rewind(10)
		
        !Load beta data
		if (ibeta==1) then !Beta(-w;w,0) case
			call loclabel(10,"-- Beta(-w,w,0) frequency",ifound)
			do i=1,ifreq
				call loclabel(10,"-- Beta(-w,w,0) frequency",ifound,0)
				read(10,*)
			end do
			read(10,*)
			read(10,*) c200tmp,beta(1,1,1)
			read(10,*) c200tmp,beta(2,1,1)
			read(10,*) c200tmp,beta(2,2,1)
			read(10,*) c200tmp,beta(3,1,1)
			read(10,*) c200tmp,beta(3,2,1)
			read(10,*) c200tmp,beta(3,3,1)
			beta(1,2,1)=beta(2,1,1)
			beta(1,3,1)=beta(3,1,1)
			beta(2,3,1)=beta(3,2,1)
			read(10,*) c200tmp,beta(1,1,2)
			read(10,*) c200tmp,beta(2,1,2)
			read(10,*) c200tmp,beta(2,2,2)
			read(10,*) c200tmp,beta(3,1,2)
			read(10,*) c200tmp,beta(3,2,2)
			read(10,*) c200tmp,beta(3,3,2)
			beta(1,2,2)=beta(2,1,2)
			beta(1,3,2)=beta(3,1,2)
			beta(2,3,2)=beta(3,2,2)
			read(10,*) c200tmp,beta(1,1,3)
			read(10,*) c200tmp,beta(2,1,3)
			read(10,*) c200tmp,beta(2,2,3)
			read(10,*) c200tmp,beta(3,1,3)
			read(10,*) c200tmp,beta(3,2,3)
			read(10,*) c200tmp,beta(3,3,3)
			beta(1,2,3)=beta(2,1,3)
			beta(1,3,3)=beta(3,1,3)
			beta(2,3,3)=beta(3,2,3)
			write(ides,"(a,/)") " Note: It is well known that the sign of hyperpolarizability of Gaussian should be multiplied by -1, the outputs below have already been corrected"
			beta=-beta
            beta_org=beta
			!Convert to other unit
			if (iunit==2) beta=beta*3.20636D-53
			if (iunit==3) beta=beta*8.63922D-33
			write(ides,*) "Frequency-dependent first hyperpolarizability Beta(-w;w,0)"
			write(ides,form) " XXX=     ",beta(1,1,1)
			write(ides,form) " XYX= YXX=",beta(1,2,1)
			write(ides,form) " YYX=     ",beta(2,2,1)
			write(ides,form) " XZX= ZXX=",beta(1,3,1)
			write(ides,form) " YZX= ZYX=",beta(2,3,1)
			write(ides,form) " ZZX=     ",beta(3,3,1)
			write(ides,form) " XXY=     ",beta(1,1,2)
			write(ides,form) " XYY= YXY=",beta(1,2,2)
			write(ides,form) " YYY=     ",beta(2,2,2)
			write(ides,form) " XZY= ZXY=",beta(1,3,2)
			write(ides,form) " YZY= ZYY=",beta(2,3,2)
			write(ides,form) " ZZY=     ",beta(3,3,2)
			write(ides,form) " XXZ=     ",beta(1,1,3)
			write(ides,form) " XYZ= YXZ=",beta(1,2,3)
			write(ides,form) " YYZ=     ",beta(2,2,3)
			write(ides,form) " XZZ= ZXZ=",beta(1,3,3)
			write(ides,form) " YZZ= ZYZ=",beta(2,3,3)
			write(ides,form) " ZZZ=     ",beta(3,3,3)
			
		else if (ibeta==2) then !used DCSHG, parsing Beta(-2w;w,w)
			call loclabel(10,"-- Beta(w,w,-2w) frequency",ifound)
            if (ifound==0) then
                write(*,"(a)") " Error: Unable to locate Beta(-2w;w,w) information. Please make sure that you have used polar=DCSHG keyword!"
                write(*,*) "Press ENTER button to exit"
                read(*,*)
                return
            end if
			do i=1,ifreq
				call loclabel(10,"-- Beta(w,w,-2w) frequency",ifound,0)
				read(10,*)
			end do
			read(10,*)
			read(10,*) c200tmp,beta(1,1,1)
			read(10,*) c200tmp,beta(1,1,2)
			beta(1,2,1)=beta(1,1,2)
			read(10,*) c200tmp,beta(1,2,2)
			read(10,*) c200tmp,beta(1,1,3)
			beta(1,3,1)=beta(1,1,3)
			read(10,*) c200tmp,beta(1,2,3)
			beta(1,3,2)=beta(1,2,3)
			read(10,*) c200tmp,beta(1,3,3)
			read(10,*) c200tmp,beta(2,1,1)
			read(10,*) c200tmp,beta(2,2,1)
			beta(2,1,2)=beta(2,2,1)
			read(10,*) c200tmp,beta(2,2,2)
			read(10,*) c200tmp,beta(2,1,3)
			beta(2,3,1)=beta(2,1,3)
			read(10,*) c200tmp,beta(2,2,3)
			beta(2,3,2)=beta(2,2,3)
			read(10,*) c200tmp,beta(2,3,3)
			read(10,*) c200tmp,beta(3,1,1)
			read(10,*) c200tmp,beta(3,1,2)
			beta(3,2,1)=beta(3,1,2)
			read(10,*) c200tmp,beta(3,2,2)
			read(10,*) c200tmp,beta(3,1,3)
			beta(3,3,1)=beta(3,1,3)
			read(10,*) c200tmp,beta(3,2,3)
			beta(3,3,2)=beta(3,2,3)
			read(10,*) c200tmp,beta(3,3,3)
			write(ides,"(a,/)") " Note: It is well known that the sign of hyperpolarizability of Gaussian should be multiplied by -1, the outputs below have already been corrected"
			beta=-beta
            beta_org=beta
			!Convert to other unit
			if (iunit==2) beta=beta*3.20636D-53
			if (iunit==3) beta=beta*8.63922D-33
			if (ibeta==2) write(ides,*) "Frequency-dependent first hyperpolarizability Beta(-2w;w,w)"
			write(ides,form) " XXX=     ",beta(1,1,1)
			write(ides,form) " YXX=     ",beta(2,1,1)
			write(ides,form) " ZXX=     ",beta(3,1,1)
			write(ides,form) " XYX= XXY=",beta(1,2,1)
			write(ides,form) " YYX= YXY=",beta(2,2,1)
			write(ides,form) " ZYX= ZXY=",beta(3,2,1)
			write(ides,form) " XYY=     ",beta(1,2,2)
			write(ides,form) " YYY=     ",beta(2,2,2)
			write(ides,form) " ZYY=     ",beta(3,2,2)
			write(ides,form) " XZX= XXZ=",beta(1,3,1)
			write(ides,form) " YZX= YXZ=",beta(2,3,1)
			write(ides,form) " ZZX= ZXZ=",beta(3,3,1)
			write(ides,form) " XZY= XYZ=",beta(1,3,2)
			write(ides,form) " YZY= YYZ=",beta(2,3,2)
			write(ides,form) " ZZY= ZYZ=",beta(3,3,2)
			write(ides,form) " XZZ=     ",beta(1,3,3)
			write(ides,form) " YZZ=     ",beta(2,3,3)
			write(ides,form) " ZZZ=     ",beta(3,3,3)
			write(ides,*)
		end if
		
        !Output beta of various forms
		betaX=0
		betaY=0
		betaZ=0
		do j=1,3
			betaX=betaX+(beta(1,j,j)+beta(j,j,1)+beta(j,1,j))/3
			betaY=betaY+(beta(2,j,j)+beta(j,j,2)+beta(j,2,j))/3
			betaZ=betaZ+(beta(3,j,j)+beta(j,j,3)+beta(j,3,j))/3
		end do
		if (iunit==1) then
			write(ides,"(' Beta_X=',f15.5,'  Beta_Y=',f15.5,'  Beta_Z=',f15.5)") betaX,betaY,betaZ
		else
			write(ides,"(' Beta_X=',1PE15.5,'  Beta_Y=',1PE15.5,'  Beta_Z=',1PE15.5)") betaX,betaY,betaZ
		end if
		write(ides,form) " Magnitude of first hyperpolarizability:",dsqrt(betaX**2+betaY**2+betaZ**2)
		betaprj=(betaX*dipole(1)+betaY*dipole(2)+betaZ*dipole(3))/dipnorm
		write(ides,form) " Projection of beta on dipole moment:",betaprj
		write(ides,form) " Beta ||     :",betaprj*3D0/5D0
		write(ides,form) " Beta ||(z)  :",betaZ*3D0/5D0
		beta_per=0
		do j=1,3
			beta_per=beta_per+(2*beta(3,j,j)+2*beta(j,j,3)-3*beta(j,3,j))/5
		end do
		write(ides,form) " Beta _|_(z) :",beta_per
        write(ides,"(/,a)") " Note: The beta information printed above corresponds to standard orientation"
        
        !Beta(-2w;w,w) as been parsed before, in this case we also output HRS related quantites
        if (ibeta==2) then
            !Calculate <beta_ZZZ^2> and <beta_XZZ^2>, without Kleinman condition approximation
            !Below, the Eqs. 4 and 5 in Phys. Chem. Chem. Phys., 10, 6223–6232 (2008) are employed
            !Note that <beta_xzz^2> is equivalent to the <beta_zxx^2> occured in many literatures
            u1=0;u2=0;u3=0;u4=0;u5=0;u6=0;u7=0;u8=0;u9=0;u10=0;u11=0
            t1=0;t2=0;t3=0;t4=0;t5=0;t6=0;t7=0;t8=0;t9=0;t10=0;t11=0
            do i=1,3
                t1=t1+beta(i,i,i)**2
                u1=u1+beta(i,i,i)**2
            end do
            do i=1,3
                do j=1,3
                    if (j/=i) then
                        t2=t2+beta(i,i,j)**2
                        t3=t3+beta(i,i,i)*beta(i,j,j)
                        t4=t4+beta(j,i,i)*beta(i,i,j)
                        t5=t5+beta(i,i,i)*beta(j,j,i)
                        t6=t6+beta(j,i,i)**2
                        u2=u2+beta(i,i,i)*beta(i,j,j)
                        u3=u3+beta(i,i,i)*beta(j,j,i)
                        u4=u4+beta(i,i,j)**2
                        u5=u5+beta(i,j,j)**2
                        u6=u6+beta(i,i,j)*beta(j,i,i)
                    end if
                    do k=1,3
                        if (k/=j.and.j/=i.and.k/=i) then
                            t7=t7+beta(i,i,j)*beta(j,k,k)
                            t8=t8+beta(j,i,i)*beta(j,k,k)
                            t9=t9+beta(i,i,j)*beta(k,k,j)
                            t10=t10+beta(i,j,k)**2
                            t11=t11+beta(i,j,k)*beta(j,i,k)
                            u7=u7+beta(i,j,j)*beta(i,k,k)
                            u8=u8+beta(i,i,k)*beta(j,j,k)
                            u9=u9+beta(i,i,j)*beta(j,k,k)
                            u10=u10+beta(i,j,k)**2
                            u11=u11+beta(i,j,k)*beta(j,i,k)
                        end if
                    end do
                end do
            end do
            betazzz2_avg = t1*(1D0/7) + t2*(4D0/35) + t3*(2D0/35) + t4*(4D0/35) + t5*(4D0/35) + &
            t6*(1D0/35) + t7*(4D0/105) + t8*(1D0/105) + t9*(4D0/105) + t10*(2D0/105) + t11*(4D0/105)
            betaxzz2_avg = u1*(1D0/35) + u2*(4D0/105) - u3*(2D0/35) + u4*(8D0/105) + u5*(3D0/35) &
            - u6*(2D0/35) + u7*(1D0/35) - u8*(2D0/105) - u9*(2D0/105) + u10*(2D0/35) - u11*(2D0/105)
        
            write(ides,*) "Note: Kleinman's symmetry condition is not employed for below quantities:"
            write(ides,"(' <beta_ZZZ^2>:',1PE16.8)") betazzz2_avg
            write(ides,"(' <beta_XZZ^2>:',1PE16.8)") betaxzz2_avg
            if (iunit==1) then
                write(ides,"(' Hyper-Rayleigh scattering (beta_HRS):',f14.3)") dsqrt(betazzz2_avg+betaxzz2_avg)
            else
                write(ides,"(' Hyper-Rayleigh scattering (beta_HRS):',1PE14.4)") dsqrt(betazzz2_avg+betaxzz2_avg)
            end if
            DRval=betazzz2_avg/betaxzz2_avg
            write(ides,"(' Depolarization ratio (DR):',f8.3)") DRval
            if (DRval<1.5D0) then
                write(ides,"(/,a,/)") " Warning: DR is smaller than 1.5, implying second-harmonic generation (SHG) resonance may occur, &
                namely half of wavelength of incident light is close to excitation energy of present system. You may consider other wavelengths of incident light"
            end if
            betaJ1norm=dsqrt(6*betazzz2_avg-9*betaxzz2_avg)
            betaJ3norm=dsqrt((-7*betazzz2_avg+63*betaxzz2_avg)/2)
            rho=betaJ3norm/betaJ1norm
            if (iunit==1) then
                write(ides,"(' |<beta J=1>|:',f14.3)") betaJ1norm
                write(ides,"(' |<beta J=3>|:',f14.3)") betaJ3norm
            else
                write(ides,"(' |<beta J=1>|:',1PE14.4)") betaJ1norm
                write(ides,"(' |<beta J=3>|:',1PE14.4)") betaJ3norm
            end if
            write(ides,"(' Nonlinear anisotropy parameter (rho):',f8.3)") rho
            write(ides,"(' Dipolar contribution to beta, phi_beta(J=1):  ',f8.3)") 1D0/(1D0+rho)
            write(ides,"(' Octupolar contribution to beta, phi_beta(J=3):',f8.3)") rho/(1D0+rho)
            write(ides,"(' < (beta_ZXZ+beta_ZZX)^2 - 2*betaZZZ*betaZXX >:',1PE16.8)") 7*betaxzz2_avg-betazzz2_avg
            !When DR<1.5: Please check whether the half of the incident light freqency (lambda/2)
            !is near to the absorption band, which is the direct reason of the resonance for SFG. (text from NLO calculator)
        
            write(*,"(/,a)") " Do you want to calculate variation of scattering intensities with respect to polarization angle of incident light? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y'.or.selectyn=='Y') then
                write(*,*) "Input initial angle in degree, e.g. -180"
                read(*,*) omegaval
                open(12,file="HRS_angle.txt",status="replace")
                do iang=1,360
                    rad=omegaval/180*pi
                    rinten=betaxzz2_avg*cos(rad)**4+betazzz2_avg*sin(rad)**4+sin(rad)**2*cos(rad)**2*(7*betaxzz2_avg-betazzz2_avg)
                    if (iunit==1) then
                        write(12,"(f6.1,f20.6)") omegaval,rinten
                    else
                        write(12,"(f6.1,1PE16.6)") omegaval,rinten
                    end if
                    omegaval=omegaval+1
                end do
                close(12)
                write(*,"(a)") " Scattering intensities corresponding to different polarization angles (degree) have been exported to HRS_angle.txt in current folder"
            end if
        end if
	end if
    
    if (ioutpol==1) then
        open(12,file="beta.txt",status="replace")
        write(12,"(1PE18.8E3)") (((beta_org(i,j,k),k=1,3),j=1,3),i=1,3)
        close(12)
        write(*,*)
        write(*,*) "Beta tensor has been exported to beta.txt in current folder"
        write(*,"(a)") " The sequence is (((beta(i,j,k),k=1,3),j=1,3),i=1,3), where k loops first. The unit is a.u."
    end if
end if


!!!! Print second polarizability part (gamma)
if (isel==7) then
    form=formother !Because magnitude of gamma is often quite large
    if (irdfreq==0.or.(irdfreq==1.and.ifreq==1)) then
        call loclabel(10,"Gamma(0;0,0,0)",ifound)
        if (ifound==0) then
            write(*,*) "Error: Unable to locate ""Gamma(0;0,0,0)""!"
            return
        end if
        read(10,*);read(10,*);read(10,*);read(10,*)
        read(10,*) c200tmp,gamma(1,1,1,1)
        read(10,*) c200tmp,gamma(1,1,2,1)
        read(10,*) c200tmp,gamma(1,1,2,2)
        read(10,*) c200tmp,gamma(2,1,2,2)
        read(10,*) c200tmp,gamma(2,2,2,2)
        read(10,*) c200tmp,gamma(1,1,3,1)
        read(10,*) c200tmp,gamma(1,1,3,2)
        read(10,*) c200tmp,gamma(2,1,3,2)
        read(10,*) c200tmp,gamma(2,2,3,2)
        read(10,*) c200tmp,gamma(1,1,3,3)
        read(10,*) c200tmp,gamma(2,1,3,3)
        read(10,*) c200tmp,gamma(2,2,3,3)
        read(10,*) c200tmp,gamma(3,1,3,3)
        read(10,*) c200tmp,gamma(3,2,3,3)
        read(10,*) c200tmp,gamma(3,3,3,3)
        do ia=1,3
            do ib=1,3
                do ic=1,3
                    do id=1,3
                        nxt=0;nyt=0;nzt=0
                        if (ia==1) nxt=nxt+1
                        if (ib==1) nxt=nxt+1
                        if (ic==1) nxt=nxt+1
                        if (id==1) nxt=nxt+1
                        if (ia==2) nyt=nyt+1
                        if (ib==2) nyt=nyt+1
                        if (ic==2) nyt=nyt+1
                        if (id==2) nyt=nyt+1
                        if (ia==3) nzt=nzt+1
                        if (ib==3) nzt=nzt+1
                        if (ic==3) nzt=nzt+1
                        if (id==3) nzt=nzt+1
                        if (nxt==3.and.nyt==1.and.nzt==0) gamma(ia,ib,ic,id)=gamma(1,1,2,1)
                        if (nxt==2.and.nyt==2.and.nzt==0) gamma(ia,ib,ic,id)=gamma(1,1,2,2)
                        if (nxt==1.and.nyt==3.and.nzt==0) gamma(ia,ib,ic,id)=gamma(2,1,2,2)
                        if (nxt==3.and.nyt==0.and.nzt==1) gamma(ia,ib,ic,id)=gamma(1,1,3,1)
                        if (nxt==2.and.nyt==1.and.nzt==1) gamma(ia,ib,ic,id)=gamma(1,1,3,2)
                        if (nxt==1.and.nyt==2.and.nzt==1) gamma(ia,ib,ic,id)=gamma(2,1,3,2)
                        if (nxt==0.and.nyt==3.and.nzt==1) gamma(ia,ib,ic,id)=gamma(2,2,3,2)
                        if (nxt==2.and.nyt==0.and.nzt==2) gamma(ia,ib,ic,id)=gamma(1,1,3,3)
                        if (nxt==1.and.nyt==1.and.nzt==2) gamma(ia,ib,ic,id)=gamma(2,1,3,3)
                        if (nxt==0.and.nyt==2.and.nzt==2) gamma(ia,ib,ic,id)=gamma(2,2,3,3)
                        if (nxt==1.and.nyt==0.and.nzt==3) gamma(ia,ib,ic,id)=gamma(3,1,3,3)
                        if (nxt==0.and.nyt==1.and.nzt==3) gamma(ia,ib,ic,id)=gamma(3,2,3,3)
                    end do
                end do
            end do
        end do
        gamma_org=gamma
		!Convert to other unit
		if (iunit==2) gamma=gamma*6.23538D-65
		if (iunit==3) gamma=gamma*5.03670D-40
        write(ides,form) " XXXX=",gamma(1,1,1,1)
        write(ides,form) " YYYY=",gamma(2,2,2,2)
        write(ides,form) " ZZZZ=",gamma(3,3,3,3)
        write(ides,form) " XXXY=XXYX=XYXX=YXXX=",gamma(1,1,2,1)
        write(ides,form) " XXYY=XYXY=XYYX=YXXY=YXYX=YYXX=",gamma(1,1,2,2)
        write(ides,form) " XYYY=YXYY=YYXY=YYYX=",gamma(2,1,2,2)
        write(ides,form) " XXXZ=XXZX=XZXX=ZXXX=",gamma(1,1,3,1)
        write(ides,form) " XXYZ=XXZY=XYXZ=XYZX=XZXY=XZYX=YXXZ=YXZX=YZXX=ZXXY=ZXYX=ZYXX=",gamma(1,1,3,2)
        write(ides,form) " XYYZ=XYZY=XZYY=YXYZ=YXZY=YYXZ=YYZX=YZXY=YZYX=ZXYY=ZYXY=ZYYX=",gamma(2,1,3,2)
        write(ides,form) " YYYZ=YYZY=YZYY=ZYYY=",gamma(2,2,3,2)
        write(ides,form) " XXZZ=XZXZ=XZZX=ZXXZ=ZXZX=ZZXX=",gamma(1,1,3,3)
        write(ides,form) " XYZZ=XZYZ=XZZY=YXZZ=YZXZ=YZZX=ZXYZ=ZXZY=ZYXZ=ZYZX=ZZXY=ZZYX=",gamma(2,1,3,3)
        write(ides,form) " YYZZ=YZYZ=YZZY=ZYYZ=ZYZY=ZZYY=",gamma(2,2,3,3)
        write(ides,form) " XZZZ=ZXZZ=ZZXZ=ZZZX=",gamma(3,1,3,3)
        write(ides,form) " YZZZ=ZYZZ=ZZYZ=ZZZY=",gamma(3,2,3,3)
        
    else if (irdfreq==1) then
        write(*,*) "Loading which type of second hyperpolarizability?"
		write(*,*) "1: Gamma(-w;w,0,0)   2: Gamma(-2w;w,w,0)"
		read(*,*) igamma
		rewind(10)
        if (igamma==1) then
			do i=1,ifreq-1
				call loclabel(10,"Gamma(-w;w,0,0)",ifound,0)
				read(10,*)
			end do
            if (ifound==0) then
                write(*,*) "Error: Unable to locate ""Gamma(-w;w,0,0)""!"
                return
            end if
            read(10,*);read(10,*);read(10,*)
            read(10,*) c200tmp,gamma(1,1,1,1)
            read(10,*) c200tmp,gamma(2,1,1,1)
            read(10,*) c200tmp,gamma(2,2,1,1)
            read(10,*) c200tmp,gamma(3,1,1,1)
            read(10,*) c200tmp,gamma(3,2,1,1)
            read(10,*) c200tmp,gamma(3,3,1,1)
            read(10,*) c200tmp,gamma(1,1,2,1)
            read(10,*) c200tmp,gamma(2,1,2,1)
            read(10,*) c200tmp,gamma(2,2,2,1)
            read(10,*) c200tmp,gamma(3,1,2,1)
            read(10,*) c200tmp,gamma(3,2,2,1)
            read(10,*) c200tmp,gamma(3,3,2,1)
            read(10,*) c200tmp,gamma(1,1,2,2)
            read(10,*) c200tmp,gamma(2,1,2,2)
            read(10,*) c200tmp,gamma(2,2,2,2)
            read(10,*) c200tmp,gamma(3,1,2,2)
            read(10,*) c200tmp,gamma(3,2,2,2)
            read(10,*) c200tmp,gamma(3,3,2,2)
            read(10,*) c200tmp,gamma(1,1,3,1)
            read(10,*) c200tmp,gamma(2,1,3,1)
            read(10,*) c200tmp,gamma(2,2,3,1)
            read(10,*) c200tmp,gamma(3,1,3,1)
            read(10,*) c200tmp,gamma(3,2,3,1)
            read(10,*) c200tmp,gamma(3,3,3,1)
            read(10,*) c200tmp,gamma(1,1,3,2)
            read(10,*) c200tmp,gamma(2,1,3,2)
            read(10,*) c200tmp,gamma(2,2,3,2)
            read(10,*) c200tmp,gamma(3,1,3,2)
            read(10,*) c200tmp,gamma(3,2,3,2)
            read(10,*) c200tmp,gamma(3,3,3,2)
            read(10,*) c200tmp,gamma(1,1,3,3)
            read(10,*) c200tmp,gamma(2,1,3,3)
            read(10,*) c200tmp,gamma(2,2,3,3)
            read(10,*) c200tmp,gamma(3,1,3,3)
            read(10,*) c200tmp,gamma(3,2,3,3)
            read(10,*) c200tmp,gamma(3,3,3,3)
            gamma(1,2,1,1)=gamma(2,1,1,1)
            gamma(1,3,1,1)=gamma(3,1,1,1)
            gamma(2,3,1,1)=gamma(3,2,1,1)
            gamma(1,1,1,2)=gamma(1,1,2,1)
            gamma(2,2,1,2)=gamma(2,2,2,1)
            gamma(3,3,1,2)=gamma(3,3,2,1)
            gamma(1,2,2,2)=gamma(2,1,2,2)
            gamma(1,3,2,2)=gamma(3,1,2,2)
            gamma(2,3,2,2)=gamma(3,2,2,2)
            gamma(1,1,1,3)=gamma(1,1,3,1)
            gamma(2,2,1,3)=gamma(2,2,3,1)
            gamma(3,3,1,3)=gamma(3,3,3,1)
            gamma(1,1,2,3)=gamma(1,1,3,2)
            gamma(2,2,2,3)=gamma(2,2,3,2)
            gamma(3,3,2,3)=gamma(3,3,3,2)
            gamma(1,2,3,3)=gamma(2,1,3,3)
            gamma(1,3,3,3)=gamma(3,1,3,3)
            gamma(2,3,3,3)=gamma(3,2,3,3)
            tmp=gamma(2,1,2,1);gamma(1,2,2,1)=tmp;gamma(2,1,1,2)=tmp;gamma(1,2,1,2)=tmp
            tmp=gamma(3,1,2,1);gamma(1,3,2,1)=tmp;gamma(3,1,1,2)=tmp;gamma(1,3,1,2)=tmp
            tmp=gamma(3,2,2,1);gamma(2,3,2,1)=tmp;gamma(3,2,1,2)=tmp;gamma(2,3,1,2)=tmp
            tmp=gamma(2,1,3,2);gamma(1,2,3,2)=tmp;gamma(2,1,2,3)=tmp;gamma(1,2,2,3)=tmp
            tmp=gamma(3,1,3,2);gamma(1,3,3,2)=tmp;gamma(3,1,2,3)=tmp;gamma(1,3,2,3)=tmp
            tmp=gamma(3,2,3,2);gamma(2,3,3,2)=tmp;gamma(3,2,2,3)=tmp;gamma(2,3,2,3)=tmp
            tmp=gamma(2,1,3,1);gamma(1,2,3,1)=tmp;gamma(2,1,1,3)=tmp;gamma(1,2,1,3)=tmp
            tmp=gamma(3,1,3,1);gamma(1,3,3,1)=tmp;gamma(3,1,1,3)=tmp;gamma(1,3,1,3)=tmp
            tmp=gamma(3,2,3,1);gamma(2,3,3,1)=tmp;gamma(3,2,1,3)=tmp;gamma(2,3,1,3)=tmp
            gamma_org=gamma
		    !Convert to other unit
		    if (iunit==2) gamma=gamma*6.23538D-65
		    if (iunit==3) gamma=gamma*5.03670D-40
            write(ides,form) " XXXX=",gamma(1,1,1,1)
            write(ides,form) " YXXX=XYXX=",gamma(2,1,1,1)
            write(ides,form) " YYXX=",gamma(2,2,1,1)
            write(ides,form) " ZXXX=XZXX=",gamma(3,1,1,1)
            write(ides,form) " ZYXX=YZXX=",gamma(3,2,1,1)
            write(ides,form) " ZZXX=",gamma(3,3,1,1)
            write(ides,form) " XXYX=XXXY=",gamma(1,1,2,1)
            write(ides,form) " YXYX=XYYX=YXXY=XYXY=",gamma(2,1,2,1)
            write(ides,form) " YYYX=YYXY=",gamma(2,2,2,1)
            write(ides,form) " ZXYX=XZYX=ZXXY=XZXY=",gamma(3,1,2,1)
            write(ides,form) " ZYYX=YZYX=ZYXY=YZXY=",gamma(3,2,2,1)
            write(ides,form) " ZZYX=ZZXY=",gamma(3,3,2,1)
            write(ides,form) " XXYY=",gamma(1,1,2,2)
            write(ides,form) " YXYY=XYYY=",gamma(2,1,2,2)
            write(ides,form) " YYYY=",gamma(2,2,2,2)
            write(ides,form) " ZXYY=XZYY=",gamma(3,1,2,2)
            write(ides,form) " ZYYY=YZYY=",gamma(3,2,2,2)
            write(ides,form) " ZZYY=",gamma(3,3,2,2)
            write(ides,form) " XXZX=XXXZ=",gamma(1,1,3,1)
            write(ides,form) " YXZX=XYZX=YXXZ=XYXZ=",gamma(2,1,3,1)
            write(ides,form) " YYZX=YYXZ=",gamma(2,2,3,1)
            write(ides,form) " ZXZX=XZZX=ZXXZ=XZXZ=",gamma(3,1,3,1)
            write(ides,form) " ZYZX=YZZX=ZYXZ=YZXZ=",gamma(3,2,3,1)
            write(ides,form) " ZZZX=ZZXZ=",gamma(3,3,3,1)
            write(ides,form) " XXZY=XXYZ=",gamma(1,1,3,2)
            write(ides,form) " YXZY=XYZY=YXYZ=XYYZ=",gamma(2,1,3,2)
            write(ides,form) " YYZY=YYYZ=",gamma(2,2,3,2)
            write(ides,form) " ZXZY=XZZY=ZXYZ=XZYZ=",gamma(3,1,3,2)
            write(ides,form) " ZYZY=YZZY=ZYYZ=YZYZ=",gamma(3,2,3,2)
            write(ides,form) " ZZZY=ZZYZ=",gamma(3,3,3,2)
            write(ides,form) " XXZZ=",gamma(1,1,3,3)
            write(ides,form) " YXZZ=XYZZ=",gamma(2,1,3,3)
            write(ides,form) " YYZZ=",gamma(2,2,3,3)
            write(ides,form) " ZXZZ=XZZZ=",gamma(3,1,3,3)
            write(ides,form) " ZYZZ=YZZZ=",gamma(3,2,3,3)
            write(ides,form) " ZZZZ=",gamma(3,3,3,3)

        else if (igamma==2) then
			do i=1,ifreq-1
				call loclabel(10,"Gamma(-2w;w,w,0)",ifound,0)
				read(10,*)
			end do
            if (ifound==0) then
                write(*,*) "Error: Unable to locate ""Gamma(-2w;w,w,0)""!"
                return
            end if
            read(10,*);read(10,*);read(10,*)
            read(10,*) c200tmp,gamma(1,1,1,1)
            read(10,*) c200tmp,gamma(2,1,1,1)
            read(10,*) c200tmp,gamma(3,1,1,1)
            read(10,*) c200tmp,gamma(1,2,1,1)
            read(10,*) c200tmp,gamma(2,2,1,1)
            read(10,*) c200tmp,gamma(3,2,1,1)
            read(10,*) c200tmp,gamma(1,2,2,1)
            read(10,*) c200tmp,gamma(2,2,2,1)
            read(10,*) c200tmp,gamma(3,2,2,1)
            read(10,*) c200tmp,gamma(1,3,1,1)
            read(10,*) c200tmp,gamma(2,3,1,1)
            read(10,*) c200tmp,gamma(3,3,1,1)
            read(10,*) c200tmp,gamma(1,3,2,1)
            read(10,*) c200tmp,gamma(2,3,2,1)
            read(10,*) c200tmp,gamma(3,3,2,1)
            read(10,*) c200tmp,gamma(1,3,3,1)
            read(10,*) c200tmp,gamma(2,3,3,1)
            read(10,*) c200tmp,gamma(3,3,3,1)
            read(10,*) c200tmp,gamma(1,1,1,2)
            read(10,*) c200tmp,gamma(2,1,1,2)
            read(10,*) c200tmp,gamma(3,1,1,2)
            read(10,*) c200tmp,gamma(1,2,1,2)
            read(10,*) c200tmp,gamma(2,2,1,2)
            read(10,*) c200tmp,gamma(3,2,1,2)
            read(10,*) c200tmp,gamma(1,2,2,2)
            read(10,*) c200tmp,gamma(2,2,2,2)
            read(10,*) c200tmp,gamma(3,2,2,2)
            read(10,*) c200tmp,gamma(1,3,1,2)
            read(10,*) c200tmp,gamma(2,3,1,2)
            read(10,*) c200tmp,gamma(3,3,1,2)
            read(10,*) c200tmp,gamma(1,3,2,2)
            read(10,*) c200tmp,gamma(2,3,2,2)
            read(10,*) c200tmp,gamma(3,3,2,2)
            read(10,*) c200tmp,gamma(1,3,3,2)
            read(10,*) c200tmp,gamma(2,3,3,2)
            read(10,*) c200tmp,gamma(3,3,3,2)
            read(10,*) c200tmp,gamma(1,1,1,3)
            read(10,*) c200tmp,gamma(2,1,1,3)
            read(10,*) c200tmp,gamma(3,1,1,3)
            read(10,*) c200tmp,gamma(1,2,1,3)
            read(10,*) c200tmp,gamma(2,2,1,3)
            read(10,*) c200tmp,gamma(3,2,1,3)
            read(10,*) c200tmp,gamma(1,2,2,3)
            read(10,*) c200tmp,gamma(2,2,2,3)
            read(10,*) c200tmp,gamma(3,2,2,3)
            read(10,*) c200tmp,gamma(1,3,1,3)
            read(10,*) c200tmp,gamma(2,3,1,3)
            read(10,*) c200tmp,gamma(3,3,1,3)
            read(10,*) c200tmp,gamma(1,3,2,3)
            read(10,*) c200tmp,gamma(2,3,2,3)
            read(10,*) c200tmp,gamma(3,3,2,3)
            read(10,*) c200tmp,gamma(1,3,3,3)
            read(10,*) c200tmp,gamma(2,3,3,3)
            read(10,*) c200tmp,gamma(3,3,3,3)
            gamma(1,1,2,1)=gamma(1,2,1,1)
            gamma(2,1,2,1)=gamma(2,2,1,1)
            gamma(3,1,2,1)=gamma(3,2,1,1)
            gamma(1,1,3,1)=gamma(1,3,1,1)
            gamma(2,1,3,1)=gamma(2,3,1,1)
            gamma(3,1,3,1)=gamma(3,3,1,1)
            gamma(1,2,3,1)=gamma(1,3,2,1)
            gamma(2,2,3,1)=gamma(2,3,2,1)
            gamma(3,2,3,1)=gamma(3,3,2,1)
            gamma(1,1,2,2)=gamma(1,2,1,2)
            gamma(2,1,2,2)=gamma(2,2,1,2)
            gamma(3,1,2,2)=gamma(3,2,1,2)
            gamma(1,1,3,2)=gamma(1,3,1,2)
            gamma(2,1,3,2)=gamma(2,3,1,2)
            gamma(3,1,3,2)=gamma(3,3,1,2)
            gamma(1,2,3,2)=gamma(1,3,2,2)
            gamma(2,2,3,2)=gamma(2,3,2,2)
            gamma(3,2,3,2)=gamma(3,3,2,2)
            gamma(1,1,2,3)=gamma(1,2,1,3)
            gamma(2,1,2,3)=gamma(2,2,1,3)
            gamma(3,1,2,3)=gamma(3,2,1,3)
            gamma(1,1,3,3)=gamma(1,3,1,3)
            gamma(2,1,3,3)=gamma(2,3,1,3)
            gamma(3,1,3,3)=gamma(3,3,1,3)
            gamma(1,2,3,3)=gamma(1,3,2,3)
            gamma(2,2,3,3)=gamma(2,3,2,3)
            gamma(3,2,3,3)=gamma(3,3,2,3)
            gamma_org=gamma
		    !Convert to other unit
		    if (iunit==2) gamma=gamma*6.23538D-65
		    if (iunit==3) gamma=gamma*5.03670D-40
            write(ides,form) " XXXX=",     gamma(1,1,1,1)
            write(ides,form) " YXXX=",     gamma(2,1,1,1)
            write(ides,form) " ZXXX=",     gamma(3,1,1,1)
            write(ides,form) " XYXX=XXYX=",gamma(1,2,1,1)
            write(ides,form) " YYXX=YXYX=",gamma(2,2,1,1)
            write(ides,form) " ZYXX=ZXYX=",gamma(3,2,1,1)
            write(ides,form) " XYYX=",     gamma(1,2,2,1)
            write(ides,form) " YYYX=",     gamma(2,2,2,1)
            write(ides,form) " ZYYX=",     gamma(3,2,2,1)
            write(ides,form) " XZXX=XXZX=",gamma(1,3,1,1)
            write(ides,form) " YZXX=YXZX=",gamma(2,3,1,1)
            write(ides,form) " ZZXX=ZXZX=",gamma(3,3,1,1)
            write(ides,form) " XZYX=XYZX=",gamma(1,3,2,1)
            write(ides,form) " YZYX=YYZX=",gamma(2,3,2,1)
            write(ides,form) " ZZYX=ZYZX=",gamma(3,3,2,1)
            write(ides,form) " XZZX=",     gamma(1,3,3,1)
            write(ides,form) " YZZX=",     gamma(2,3,3,1)
            write(ides,form) " ZZZX=",     gamma(3,3,3,1)
            write(ides,form) " XXXY=",     gamma(1,1,1,2)
            write(ides,form) " YXXY=",     gamma(2,1,1,2)
            write(ides,form) " ZXXY=",     gamma(3,1,1,2)
            write(ides,form) " XYXY=XXYY=",gamma(1,2,1,2)
            write(ides,form) " YYXY=YXYY=",gamma(2,2,1,2)
            write(ides,form) " ZYXY=ZXYY=",gamma(3,2,1,2)
            write(ides,form) " XYYY=",     gamma(1,2,2,2)
            write(ides,form) " YYYY=",     gamma(2,2,2,2)
            write(ides,form) " ZYYY=",     gamma(3,2,2,2)
            write(ides,form) " XZXY=XXZY=",gamma(1,3,1,2)
            write(ides,form) " YZXY=YXZY=",gamma(2,3,1,2)
            write(ides,form) " ZZXY=ZXZY=",gamma(3,3,1,2)
            write(ides,form) " XZYY=XYZY=",gamma(1,3,2,2)
            write(ides,form) " YZYY=YYZY=",gamma(2,3,2,2)
            write(ides,form) " ZZYY=ZYZY=",gamma(3,3,2,2)
            write(ides,form) " XZZY=",     gamma(1,3,3,2)
            write(ides,form) " YZZY=",     gamma(2,3,3,2)
            write(ides,form) " ZZZY=",     gamma(3,3,3,2)
            write(ides,form) " XXXZ=",     gamma(1,1,1,3)
            write(ides,form) " YXXZ=",     gamma(2,1,1,3)
            write(ides,form) " ZXXZ=",     gamma(3,1,1,3)
            write(ides,form) " XYXZ=XXYZ=",gamma(1,2,1,3)
            write(ides,form) " YYXZ=YXYZ=",gamma(2,2,1,3)
            write(ides,form) " ZYXZ=ZXYZ=",gamma(3,2,1,3)
            write(ides,form) " XYYZ=",     gamma(1,2,2,3)
            write(ides,form) " YYYZ=",     gamma(2,2,2,3)
            write(ides,form) " ZYYZ=",     gamma(3,2,2,3)
            write(ides,form) " XZXZ=XXZZ=",gamma(1,3,1,3)
            write(ides,form) " YZXZ=YXZZ=",gamma(2,3,1,3)
            write(ides,form) " ZZXZ=ZXZZ=",gamma(3,3,1,3)
            write(ides,form) " XZYZ=XYZZ=",gamma(1,3,2,3)
            write(ides,form) " YZYZ=YYZZ=",gamma(2,3,2,3)
            write(ides,form) " ZZYZ=ZYZZ=",gamma(3,3,2,3)
            write(ides,form) " XZZZ=",     gamma(1,3,3,3)
            write(ides,form) " YZZZ=",     gamma(2,3,3,3)
            write(ides,form) " ZZZZ=",     gamma(3,3,3,3)
        end if
    end if
    !Print all gamma components one by one
    !write(*,*) "Debug"
    !do ia=1,3
    !    do ib=1,3
    !        do ic=1,3
    !            do id=1,3
    !                write(ides,form) ' '//lb(ia)//lb(ib)//lb(ic)//lb(id)//'=',gamma(ia,ib,ic,id)
    !            end do
    !        end do
    !    end do
    !end do
	gammaX=0;gammaY=0;gammaZ=0
	do i=1,3
		gammaX=gammaX+gamma(1,i,i,1)+gamma(1,i,1,i)+gamma(1,1,i,i)
		gammaY=gammaY+gamma(2,i,i,2)+gamma(2,i,2,i)+gamma(2,2,i,i)
		gammaZ=gammaZ+gamma(3,i,i,3)+gamma(3,i,3,i)+gamma(3,3,i,i)
	end do
	gammaX=gammaX/15;gammaY=gammaY/15;gammaZ=gammaZ/15
	gammatot=dsqrt(gammaX**2+gammaY**2+gammaZ**2)
	gammaavg1=gammaX+gammaY+gammaZ
	gammaavg2=( gamma(1,1,1,1)+gamma(2,2,2,2)+gamma(3,3,3,3) + gamma(1,1,2,2)+gamma(1,1,3,3)+gamma(2,2,3,3) &
    + gamma(2,2,1,1)+gamma(3,3,1,1)+gamma(3,3,2,2) )/5
    gamma_nor=0
    do i=1,3
        do j=1,3
            gamma_nor=gamma_nor+ 2*gamma(i,j,j,i)-gamma(i,i,j,j)
        end do
    end do
    gamma_nor=gamma_nor/15
    write(ides,*)
    write(ides,form) " Magnitude of gamma:  ",gammatot
    write(ides,form) " X component of gamma:",gammaX
    write(ides,form) " Y component of gamma:",gammaY
    write(ides,form) " Z component of gamma:",gammaZ
    write(ides,form) " Average of gamma (definition 1), gamma ||:",gammaavg1
    write(ides,form) " Average of gamma (definition 2):",gammaavg2
    write(ides,form) " gamma _|_:",gamma_nor
    write(ides,"(/,a)") " Note: Gamma (in input orientation) has been given above. If you want to obtain beta information, please use option 1"
        
    if (ioutpol==1) then
        open(12,file="gamma.txt",status="replace")
        write(12,"(1PE18.8E3)") ((((gamma_org(i,j,k,l),l=1,3),k=1,3),j=1,3),i=1,3)
        close(12)
        write(*,*)
        write(*,*) "Gamma tensor has been exported to gamma.txt in current folder"
        write(*,"(a)") " The sequence is ((((gamma(i,j,k,l),l=1,3),k=1,3),j=1,3),i=1,3), where l loops first. The unit is a.u."
    end if
end if !End gamma

close(10)
if (ides==11) close(ides)
goto 10
end subroutine





!!-----------------------------------------------------------------------------------------------------------
!!--------- Sum-over-states (SOS) calculation for (hyper)polarizability and two/three-level analyses --------
!!-----------------------------------------------------------------------------------------------------------
!Programmed based on the formulae in Sasagane et al. J. Chem. Phys., 99, 3738 (1993)
subroutine SOS
use defvar
use util
implicit real*8 (a-h,o-z)
character transmodestr*80,c80tmp*80,c200tmp*80
character :: dirlab(3)=(/ "X","Y","Z" /)
real*8,allocatable :: trandip(:,:,:) !Transition dipole moment between i and j in X,Y,Z. 0 corresponds to ground state
real*8,allocatable :: excene(:) !Excitation energy
real*8 alpha(3,3),alphatmp(3,3),beta(3,3,3),gamma(3,3,3,3),gamma1(3,3,3,3),gamma2(3,3,3,3),delta(3,3,3,3,3)
real*8 eigval(3),eigvecmat(3,3),tmpw(5)
real*8,allocatable :: freqlist(:,:) !Store the frequency to be calculated for beta and gamma
real*8 term(3)
integer tmpdir(5),arrb(6,3),arrg(24,4),arrd(120,5),dir1,dir2,dir3,dir4,dir5

write(*,*) "Loading data..."
open(10,file=filename,status="old")
call loclabel(10,"Gaussian, Inc",igauout,maxline=100)
rewind(10)
if (igauout==1) then !Load excitation energies and <0|r|n>
	write(*,*) "This is a Gaussian output file"
	call loclabel(10,"Excitation energies and oscillator strengths:")
	read(10,*)
	nstates=0 !The number of transitions
	do while(.true.)
		call loclabel(10,"Excited State",ifound,0)
		if (ifound==1) then
			nstates=nstates+1
			read(10,*)
		else
			exit
		end if
	end do
	allocate(trandip(0:nstates,0:nstates,3),excene(0:nstates))
	trandip=0
	excene=0
	call loclabel(10,"Ground to excited state transition electric dipole moments")
	read(10,*)
	read(10,*)
	do istat=1,nstates
		read(10,*) inouse,trandip(0,istat,:) !Transition dipole moment from ground state (0) to excited "istate"
	end do
	call loclabel(10,"Excitation Energies [eV] at current iteration:",ifound)
	if (ifound==1) then !Read the excitation energies shown in the last step of iteration process will be more accurate, but #P must be used
		ncyc=1
		do while(.true.) !Find the position of the last iteration
			call loclabel(10,"Excitation Energies [eV] at current iteration:",ifound,0)
			if (ifound==0) exit
			read(10,*)
			ncyc=ncyc+1
		end do
		rewind(10)
		do icyc=1,ncyc
			call loclabel(10,"Excitation Energies [eV] at current iteration:",ifound,0)
		end do
		read(10,*)
		do istat=1,nstates
			read(10,"(14x)",advance="no")
			read(10,*) excene(istat)
		end do
	else !Read excitation energies from normal output
		call loclabel(10,"Excitation energies and oscillator strengths:")
		do istat=1,nstates
			call loclabel(10,"Excited State",ifound,0)
			read(10,"(a)") transmodestr
			do i=10,70
				if (transmodestr(i:i+1)=="eV") exit
			end do
			read(transmodestr(i-10:i-1),*) excene(istat) !Read as eV
		end do
	end if
	call loclabel(10,"Dipole moment (field-independent basis, Debye)")
	read(10,*)
	read(10,*) c80tmp,xtmp,c80tmp,ytmp,c80tmp,ztmp
	trandip(0,0,1)=xtmp/au2debye
	trandip(0,0,2)=ytmp/au2debye
	trandip(0,0,3)=ztmp/au2debye
	ionlyalpha=1
else !Load excitation energies and all <m|r|n> from plain text file
	ionlyalpha=0
	read(10,*) nstates
	if (nstates<0) ionlyalpha=1 !Only calculate polarizability, not hyperpolarizability, so will not read transition dipole moments between excited states
	nstates=abs(nstates)
	allocate(trandip(0:nstates,0:nstates,3),excene(0:nstates))
	do i=1,nstates !Read as eV
		read(10,*) inouse,excene(i)
	end do
	do i=0,nstates !i-i corresponds to dipole moment of corresponding state
		do j=i,nstates
			read(10,*) inouse,inouse,trandip(i,j,:)
		end do
		if (ionlyalpha==1) exit
	end do
end if
close(10)
excene(0)=0
excene=excene/au2eV
!Transition dipole moments are loaded as upper trigonal matrix, now we convert it as symmetry matrix
do i=0,nstates
	do j=i+1,nstates
		trandip(j,i,:)=trandip(i,j,:)
	end do
end do

!! Output some information
! write(*,*) "  State#      Exc.ene.(a.u.)     Transition dipole moment in X,Y,Z (a.u.)"
! do istat=1,nstates
! 	write(*,"(i8,f20.6,3f15.6)") istat,excene(istat),trandip(0,istat,:)
! end do
do istat=1,nstates
	write(*,"(' State',i6,'   Excitation energy:',f12.5,' a.u.',f14.6,' eV')") istat,excene(istat),excene(istat)*au2eV
end do
write(*,"(' There are',i6,' excited states')") nstates
write(*,"(' Dipole moment of ground state contributed by electrons:',/,'   X=',f12.5,'   Y=',f12.5,'   Z=',f12.5,' a.u.')") trandip(0,0,:)
write(*,*)
write(*,*) "NOTE: Unless otherwise specified, all units used in this function are a.u."
!Gaussian output file is impossible to provide <m|r|n>, even if alltransitiondensities is used for CIS, it doesn't output <m|r|m>
do while(.true.)
write(*,*)
write(*,*) "  ------- Sum-over-states (SOS) calculation for (hyper)polarizability -------"
write(*,*) "0 Return"
write(*,*) "1 Calculate polarizability (alpha)"
if (ionlyalpha==0) write(*,*) "2 Calculate first hyperpolarizability (beta)"
if (ionlyalpha==0) write(*,*) "3 Calculate second hyperpolarizability (gamma)"
if (ionlyalpha==0) write(*,*) "4 Calculate third hyperpolarizability (delta)"
write(*,*) "5 Show the variation of alpha w.r.t. the number of states in consideration"
if (ionlyalpha==0) write(*,*) "6 Show the variation of beta w.r.t. the number of states in consideration"
if (ionlyalpha==0) write(*,*) "7 Show the variation of gamma w.r.t. the number of states in consideration"
write(*,*) "15 Calculate alpha in a range of frequencies"
if (ionlyalpha==0) write(*,*) "16 Calculate beta in a range of frequencies"
if (ionlyalpha==0) write(*,*) "17 Calculate gamma in a range of frequencies"
if (ionlyalpha==0) write(*,*) "19 Scanning w1 and w2 of beta(-(w1+w2);w1,w2)"
if (ionlyalpha==0) write(*,*) "20 Two or three-level model analysis of beta"
read(*,*) isel

if (isel==0) then
	return
    
else if (isel==1.or.isel==5.or.isel==15) then !Analysis of polarizability (alpha)
!1=Calculate polarizability
!5=Show the variation of alpha w.r.t. the number of states in consideration
!15=Calculate alpha in a range of frequencies
	if (isel==1.or.isel==5) then
		write(*,*) "Input frequency of external field w for alpha(-w;w) in a.u., e.g. 0.25"
        write(*,*) "You can also input the value in nm, e.g. 532 nm"
		write(*,*) "Note: 0 corresponds to static case"
        read(*,"(a)") c80tmp
		read(c80tmp,*) freq
		if (index(c80tmp,"nm")/=0) freq=1239.842D0/au2eV/freq !Convert to a.u.
		if (freq/=0) then
			wavlen=1239.842D0/(freq*au2eV)
			write(*,"(' Wavelength of w:',f12.6,' a.u.',f12.3,' nm',/)") freq,wavlen
        else
           write(*,*) "External field is static" 
		end if
	else if (isel==15) then
		write(*,*) "Input lower, upper limits and stepsize of w for alpha(-w;w)"
		write(*,*) "e.g. 0.2,0.5,0.01"
		read(*,*) freq,freqend,freqstep
	end if
	
	if (isel==1) then
		istart=nstates
		iend=nstates
	else if (isel==5) then
		istart=1
		iend=nstates
		open(10,file="alpha_n.txt",status="replace")
	else if (isel==15) then
		istart=nstates
		iend=nstates
		open(10,file="alpha_w.txt",status="replace")
	end if
	
	write(*,*) "Please wait..."
	do numstat=istart,iend !Cycle number of states
		do while(.true.) !For isel==15, vary frequency (freq) from initial value until reaching ending (freqend); For other cases, do once
			do idir=1,3
				do jdir=1,3
					tmpval=0
					do istat=1,numstat
						tmpval1=trandip(0,istat,idir)*trandip(istat,0,jdir)/(excene(istat)-freq)
						tmpval2=trandip(0,istat,jdir)*trandip(istat,0,idir)/(excene(istat)+freq)
						tmpval=tmpval+tmpval1+tmpval2
					end do
					alpha(idir,jdir)=tmpval
				end do
			end do
			alphaiso=(alpha(1,1)+alpha(2,2)+alpha(3,3))/3D0
			term1=(alpha(1,1)-alpha(2,2))**2 + (alpha(1,1)-alpha(3,3))**2 + (alpha(2,2)-alpha(3,3))**2
			term2=6*(alpha(1,2)**2+alpha(1,3)**2+alpha(2,3)**2)
			alphaani1=dsqrt((term1+term2)/2D0)
            alphatmp=alpha
			call diagmat(alphatmp,eigvecmat,eigval,300,1D-10)
			call sort(eigval)
			alphaani2=eigval(3)-(eigval(1)+eigval(2))/2D0
		
			if (isel==1) then
				write(*,*) "Polarizability tensor:"
				write(*,*) "             1              2              3"
				do idir=1,3
					write(*,"(i3,3f15.6)") idir,alpha(idir,:)
				end do
				write(*,"(' Isotropic average polarizability:',f15.6)") alphaiso
				write(*,"(' Isotropic average polarizability volume:',f15.6,' Angstrom^3')") alphaiso*0.14818470D0
				write(*,"(' Polarizability anisotropy (definition 1):',f15.6)") alphaani1
				write(*,"(' Eigenvalues:',3f15.6)") eigval(:)
				write(*,"(' Polarizability anisotropy (definition 2):',f15.6)") alphaani2
			else if (isel==5) then
				write(10,"(i6,9f15.6)") numstat,alphaiso,alphaani1,alphaani2,alpha(1,1),alpha(2,2),alpha(3,3),alpha(1,2),alpha(1,3),alpha(2,3)
			else if (isel==15) then
				write(10,"(f12.6,9(1PE14.5))") freq,alphaiso,alphaani1,alphaani2,alpha(1,1),alpha(2,2),alpha(3,3),alpha(1,2),alpha(1,3),alpha(2,3)
			end if
            if (isel==1.or.isel==5) then
                exit
            else if (isel==15) then
			    freq=freq+freqstep
			    if (freq>freqend) exit
            end if
		end do
	end do
	if (isel==5.or.isel==15) then
		close(10)
		if (isel==5) write(*,*) "Done! The result has been outputted to alpha_n.txt in current folder"
		if (isel==15) write(*,*) "Done! The result has been outputted to alpha_w.txt in current folder"
		write(*,*) "The correspondence between columns and information in this file is as follows"
		if (isel==5) write(*,*) "Column 1:  The number of states in consideration"
		if (isel==15)  write(*,*) "Column 1:  Frequency of external field"
		write(*,*) "Column 2:  Isotropic average polarizability"
		write(*,*) "Column 3:  Polarizability anisotropy (definition 1)"
		write(*,*) "Column 4:  Polarizability anisotropy (definition 2)"
		write(*,*) "Column 5:  XX"
		write(*,*) "Column 6:  YY"
		write(*,*) "Column 7:  ZZ"
		write(*,*) "Column 8:  XY"
		write(*,*) "Column 9:  XZ"
		write(*,*) "Column 10: YZ"
	end if

else if (isel==2.or.isel==6.or.isel==16.or.isel==19) then !Analysis of first hyperpolarizability (beta)
!2=Calculate first hyperpolarizability (beta)
!6=Show the variation of beta w.r.t. the number of states in consideration
!16=Calculate beta in a range of frequencies
!19=Scanning w1 and w2 of beta(-(w1+w2);w1,w2)
    !Load frequency setting
	if (allocated(freqlist)) deallocate(freqlist)
	if (isel==2.or.isel==6) then
		nfreq=1 !Only one frequency is considered
		allocate(freqlist(1,2))
		write(*,*) "Input w1 and w2 in a.u. for beta(-(w1+w2);w1,w2), e.g. 0.25,-0.13"
        write(*,*) "You can also input the values in nm, e.g. 450,-532 nm"
        read(*,"(a)") c80tmp
		read(c80tmp,*) freqlist(1,:)
		if (index(c80tmp,"nm")/=0) freqlist=1239.842D0/au2eV/freqlist
		if (freqlist(1,1)/=0) then
			wavlen1=1239.842D0/(freqlist(1,1)*au2eV)
			write(*,"(' Wavelength of w1:',f12.6,' a.u.',f12.3,' nm')") freqlist(1,1),wavlen1
        else
            write(*,*) "w1 is static"
		end if
		if (freqlist(1,2)/=0) then
			wavlen2=1239.842D0/(freqlist(1,2)*au2eV)
			write(*,"(' Wavelength of w2:',f12.6,' a.u.',f12.3,' nm')") freqlist(1,2),wavlen2
        else
            write(*,*) "w2 is static"
		end if
		write(*,*)
	else if (isel==16) then
		write(*,*) "Input the file recording frequency list, e.g. C:\freqlist.txt"
		write(*,"(a)") " The file should contain two columns, corresponding to frequency of w1 and w2 in a.u., respectively"
		do while(.true.)
			read(*,"(a)") c200tmp
			inquire(file=c200tmp,exist=alive)
			if (alive) exit
			write(*,*) "Cannot find the file, input again"
		end do
		open(10,file=c200tmp,status="old")
		nfreq=totlinenum(10,1)
		allocate(freqlist(nfreq,2))
		do ifreq=1,nfreq
			read(10,*) freqlist(ifreq,:)
		end do
		close(10)
		write(*,*) "The frequencies loaded:"
		do ifreq=1,nfreq
			write(*,"(' #',i5,'  w1=',f10.5,' a.u.',f10.3,' nm   w2=',f10.5' a.u.',f10.3,' nm')") &
			ifreq,freqlist(ifreq,1),1239.842D0/(freqlist(ifreq,1)*au2eV),freqlist(ifreq,2),1239.842D0/(freqlist(ifreq,2)*au2eV)
		end do
		write(*,*)
	else if (isel==19) then
        write(*,*) "Input initial value, ending value and number of steps for w1 (in a.u.)"
        write(*,*) "Example: -0.15,0.15,10"
        read(*,*) begw1,endw1,nw1
        stepw1=(endw1-begw1)/nw1
        write(*,"(' Stepsize of w1:',f12.6,' a.u.')") stepw1
        write(*,*) "Input initial value, ending value and number of steps for w2 (in a.u.)"
        write(*,*) "Example: -0.1,0.1,10"
        read(*,*) begw2,endw2,nw2
        stepw2=(endw2-begw2)/nw2
        write(*,"(' Stepsize of w2:',f12.6,' a.u.')") stepw2
        nfreq=(nw1+1)*(nw2+1)
		allocate(freqlist(nfreq,2))
        ifreq=0
        do iw1=0,nw1
            do iw2=0,nw2
                ifreq=ifreq+1
                freqlist(ifreq,1)=begw1+iw1*stepw1
                freqlist(ifreq,2)=begw2+iw2*stepw2
            end do
        end do
	end if
    
    !Initialize variables
	if (isel==2) then
		istart=nstates
		iend=nstates
	else if (isel==6) then
		istart=1
		iend=nstates
		open(10,file="beta_n.txt",status="replace")
		open(11,file="beta_n_comp.txt",status="replace")
        !Write first line of beta_n_comp.txt, namely meaning of each column
        write(11,"(a)",advance="no") "  N_states "
		do idir=1,3
			do jdir=1,3
				do kdir=1,3
                    write(11,"(5x,a,6x)",advance="no") dirlab(idir)//dirlab(jdir)//dirlab(kdir)
                end do
            end do
        end do
        write(11,*)
	else if (isel==16.or.isel==19) then
		istart=nstates
		iend=nstates
		open(10,file="beta_w.txt",status="replace")
		open(11,file="beta_w_comp.txt",status="replace")
        !Write first line of beta_w_comp.txt, namely meaning of each column
        write(11,"(a)",advance="no") "    Freq. 1     Freq. 2  "
		do idir=1,3
			do jdir=1,3
				do kdir=1,3
                    write(11,"(5x,a,6x)",advance="no") dirlab(idir)//dirlab(jdir)//dirlab(kdir)
                end do
            end do
        end do
        write(11,*)
	end if
	
	write(*,*) "Please wait..."
    if (isel==16.or.isel==19) call showprog(0,nfreq)
	call fullarrange(arrb,6,3) !Generate full arrangement matrix (3!=6 :3) for beta, each row corresponds to one permutation, e.g. 231
	do numstat=istart,iend
	do ifreq=1,nfreq
        if (isel==16.or.isel==19) call showprog(ifreq,nfreq)
		freq1=freqlist(ifreq,1)
		freq2=freqlist(ifreq,2)
		freqtot=freq1+freq2
		tmpw(1:3)=(/ -freqtot,freq1,freq2 /)
		do idir=1,3
			do jdir=1,3
				do kdir=1,3

					tmpval=0
					tmpdir(1:3)=(/ idir,jdir,kdir /)
					do iper=1,6 !Do permutation, arrb(1,:)=1,2,3
						dir1=tmpdir(arrb(iper,1))
						dir2=tmpdir(arrb(iper,2))
						dir3=tmpdir(arrb(iper,3))
						w0=tmpw(arrb(iper,1))
						w2=tmpw(arrb(iper,3))
						do istat=1,numstat
							do jstat=1,numstat
								cen=trandip(istat,jstat,dir2)
								if (istat==jstat) cen=cen-trandip(0,0,dir2)
								tmpval=tmpval+trandip(0,istat,dir1)*cen*trandip(jstat,0,dir3)/(excene(istat)+w0)/(excene(jstat)-w2)
							end do
						end do
					end do
					beta(idir,jdir,kdir)=tmpval
					
					!Below is the code manually considering each permutation, for teaching purpose, but foolish
! 					tmpval=0
! 					do istat=1,numstat
! 						do jstat=1,numstat
! 							! Consider six permutations, i=-freqtot, j=freq1, k=freq2,  the denominator is (+1th),(-3th)
! 							! 1_3
! 							! ijk; -freqtot,-freq2
! 							! ikj; -freqtot,-freq1
! 							! jik; +freq1  ,-freq2
! 							! jki; +freq1  ,+freqtot
! 							! kij; +freq2  ,-freq1
! 							! kji; +freq2  ,+freqtot
! 							t1c=trandip(istat,jstat,jdir)
! 							t2c=trandip(istat,jstat,kdir)
! 							t3c=trandip(istat,jstat,idir)
! 							t4c=trandip(istat,jstat,kdir)
! 							t5c=trandip(istat,jstat,idir)
! 							t6c=trandip(istat,jstat,jdir)
! 							if (istat==jstat) then
! 								t1c=t1c-trandip(0,0,jdir)
! 								t2c=t2c-trandip(0,0,kdir)
! 								t3c=t3c-trandip(0,0,idir)
! 								t4c=t4c-trandip(0,0,kdir)
! 								t5c=t5c-trandip(0,0,idir)
! 								t6c=t6c-trandip(0,0,jdir)
! 							end if
! 							t1=trandip(0,istat,idir)*t1c*trandip(jstat,0,kdir) /(excene(istat)-freqtot)/(excene(jstat)-freq2)
! 							t2=trandip(0,istat,idir)*t2c*trandip(jstat,0,jdir) /(excene(istat)-freqtot)/(excene(jstat)-freq1)
! 							t3=trandip(0,istat,jdir)*t3c*trandip(jstat,0,kdir) /(excene(istat)+freq1)  /(excene(jstat)-freq2)
! 							t4=trandip(0,istat,jdir)*t4c*trandip(jstat,0,idir) /(excene(istat)+freq1)  /(excene(jstat)+freqtot)
! 							t5=trandip(0,istat,kdir)*t5c*trandip(jstat,0,jdir) /(excene(istat)+freq2)  /(excene(jstat)-freq1)
! 							t6=trandip(0,istat,kdir)*t6c*trandip(jstat,0,idir) /(excene(istat)+freq2)  /(excene(jstat)+freqtot)
! 							tmpval=tmpval+t1+t2+t3+t4+t5+t6
! 						end do
! 					end do
! 					beta(idir,jdir,kdir)=tmpval
				end do
			end do
		end do

		betaX=0
		betaY=0
		betaZ=0
		do j=1,3
			betaX=betaX+(beta(1,j,j)+beta(j,j,1)+beta(j,1,j))/3
			betaY=betaY+(beta(2,j,j)+beta(j,j,2)+beta(j,2,j))/3
			betaZ=betaZ+(beta(3,j,j)+beta(j,j,3)+beta(j,3,j))/3
		end do
		betatot=dsqrt(betaX**2+betaY**2+betaZ**2)
		dipx=trandip(0,0,1)
		dipy=trandip(0,0,2)
		dipz=trandip(0,0,3)
		dipnorm=dsqrt(dipx**2+dipy**2+dipz**2)
		betaprj=(betaX*dipx+betaY*dipy+betaZ*dipz)/dipnorm
		beta_per=0
		do j=1,3
			beta_per=beta_per+(2*beta(3,j,j)+2*beta(j,j,3)-3*beta(j,3,j))/5
		end do
		
		if (isel==2) then
			write(*,*) "First hyperpolarizability tensor:"
			do jdir=1,3
				do kdir=1,3
					do idir=1,3
						write(*,"(2x,3a,'=',f17.5,2x)",advance="no") dirlab(idir),dirlab(jdir),dirlab(kdir),beta(idir,jdir,kdir)
						if (idir==3) write(*,*)
					end do
				end do
			end do
			write(*,"(/,' Beta_X:',f17.5,'  Beta_Y:',f17.5,'  Beta_Z:',f17.5)") betaX,betaY,betaZ
			write(*,"(a,f17.5)") " Magnitude of beta:",betatot
			write(*,"(a,f17.5)") " Projection of beta on dipole moment:",betaprj
			write(*,"(a,f17.5)") " Beta ||     :",betaprj*3D0/5D0
			write(*,"(a,f17.5)") " Beta ||(z)  :",betaZ*3D0/5D0
			write(*,"(a,f17.5)") " Beta _|_(z) :",beta_per
		else if (isel==6) then
			write(10,"(i6,8(1PE14.5))") numstat,betaX,betaY,betaZ,betatot,betaprj,betaprj*3D0/5D0,betaZ*3D0/5D0,beta_per
            write(11,"(i10)",advance="no") numstat
            write(11,"(27(1PE14.5))",advance="no") (((beta(idir,jdir,kdir),kdir=1,3),jdir=1,3),idir=1,3)
			write(11,*)
		else if (isel==16.or.isel==19) then
			write(10,"(2f12.6,8(1PE14.5))") freq1,freq2,betaX,betaY,betaZ,betatot,betaprj,betaprj*3D0/5D0,betaZ*3D0/5D0,beta_per
            write(11,"(2f12.6)",advance="no") freq1,freq2
            write(11,"(27(1PE14.5))",advance="no") (((beta(idir,jdir,kdir),kdir=1,3),jdir=1,3),idir=1,3)
			write(11,*)
		end if
	end do
	end do
	
	if (isel==6.or.isel==16.or.isel==19) then
		close(10)
        close(11)
		if (isel==6) then
			write(*,*) "Done! The result has been outputted to beta_n.txt in current folder"
			write(*,*) "The correspondence between columns and information in this file is as follows"
			write(*,*) "Column 1:  The number of states in consideration"
			write(*,*) "Column 2:  Beta_X"
			write(*,*) "Column 3:  Beta_Y"
			write(*,*) "Column 4:  Beta_Z"
			write(*,*) "Column 5:  Magnitude of hyperpolarizability"
			write(*,*) "Column 6:  Hyperpolarizability component along dipole moment direction"
			write(*,*) "Column 7:  Beta ||"
			write(*,*) "Column 8:  Beta ||(z)"
			write(*,*) "Column 9:  Beta _|_(z)"
            write(*,"(a)") " In addition, all components of beta with respect to the number of considered states &
            have been exported to beta_n_comp.txt in current folder"
		else if (isel==16.or.isel==19) then
			write(*,*) "Done! The result has been outputted to beta_w.txt in current folder"
			write(*,*) "The correspondence between columns and information in this file is as follows"
			write(*,*) "Column 1:  Frequency of the first external field (w1)"
			write(*,*) "Column 2:  Frequency of the second external field (w2)"
			write(*,*) "Column 3:  Beta_X"
			write(*,*) "Column 4:  Beta_Y"
			write(*,*) "Column 5:  Beta_Z"
			write(*,*) "Column 6:  Magnitude of hyperpolarizability"
			write(*,*) "Column 7:  Hyperpolarizability component along dipole moment direction"
			write(*,*) "Column 8:  Beta ||"
			write(*,*) "Column 9:  Beta ||(z)"
			write(*,*) "Column 10: Beta _|_(z)"
            write(*,"(a)") " In addition, all components of beta with respect to frequencies &
            have been exported to beta_w_comp.txt in current folder"
		end if 
	end if
	

else if (isel==3.or.isel==7.or.isel==17) then !Analysis of second hyperpolarizability (gamma)
!3=Calculate second hyperpolarizability (gamma)
!7=Show the variation of gamma w.r.t. the number of states in consideration
!17=Calculate gamma in a range of frequencies
    !Load frequency setting
	if (allocated(freqlist)) deallocate(freqlist)
	if (isel==3.or.isel==7) then
		nfreq=1
		allocate(freqlist(1,3))
		write(*,*) "Input w1, w2, w3 in a.u. for gamma(-(w1+w2+w3);w1,w2,w3), e.g. 0.13,-0.13,0"
        write(*,*) "You can also input the values in nm, e.g. 450.6,-532,532 nm"
        read(*,"(a)") c80tmp
		read(c80tmp,*) freqlist(1,:)
		if (index(c80tmp,"nm")/=0) freqlist=1239.842D0/au2eV/freqlist
		do i=1,3
			if (freqlist(1,i)/=0) then
				wavlen=1239.842D0/(freqlist(1,i)*au2eV)
				write(*,"(' Wavelength of w',i1,':',f12.6,' a.u.',f12.3,' nm')") i,freqlist(1,i),wavlen
            else
                write(*,"(' w',i1,' is static')") i
			end if
		end do
		write(*,*)
	else if (isel==17) then
		write(*,*) "Input the file recording frequency list, e.g. C:\freqlist.txt"
		write(*,"(a)") " The file should contain three columns, corresponding to frequency of w1, w2 and w3 in a.u., respectively"
		do while(.true.)
			read(*,"(a)") c200tmp
			inquire(file=c200tmp,exist=alive)
			if (alive) exit
			write(*,*) "Cannot find the file, input again"
		end do
		open(10,file=c200tmp,status="old")
		nfreq=totlinenum(10,1)
		allocate(freqlist(nfreq,3))
		do ifreq=1,nfreq
			read(10,*) freqlist(ifreq,:)
		end do
		close(10)
		write(*,*) "The frequencies loaded:"
		do ifreq=1,nfreq
			write(*,"(' #',i5,'   w1=',f10.5,' a.u.    w2=',f10.5,' a.u.    w3=',f10.5,' a.u.')") ifreq,freqlist(ifreq,:)
		end do
		write(*,*)
	end if
	
	if (isel==3.or.isel==17) then
		write(*,"(' Consider how many states? Should be <=',i6)") nstates
		read(*,*) istart
		iend=istart
		nstatstep=1
	else if (isel==7) then
		write(*,*) "Input upper limit and stepsize of the number of considered states, e.g. 150,2"
		read(*,*) iend,nstatstep
		istart=1
		if (iend>nstates) iend=nstates
	end if
	if (isel==7) then
        open(10,file="gamma_n.txt",status="replace")
        open(11,file="gamma_n_comp.txt",status="replace")
        open(12,file="gamma_I_n_comp.txt",status="replace")
        open(13,file="gamma_II_n_comp.txt",status="replace")
        !Write first line of gamma_n_comp.txt, namely meaning of each column
        do ifile=11,13
            write(ifile,"(a)",advance="no") "  N_states "
		    do idir=1,3
			    do jdir=1,3
				    do kdir=1,3
				        do ldir=1,3
                            write(ifile,"(5x,a,5x)",advance="no") dirlab(idir)//dirlab(jdir)//dirlab(kdir)//dirlab(ldir)
                        end do
                    end do
                end do
            end do
            write(ifile,*)
        end do
	else if (isel==17) then
        open(10,file="gamma_w.txt",status="replace")
        open(11,file="gamma_w_comp.txt",status="replace")
        open(12,file="gamma_I_w_comp.txt",status="replace")
        open(13,file="gamma_II_w_comp.txt",status="replace")
        !Write first line of gamma_w_comp.txt, namely meaning of each column
        do ifile=11,13
            write(ifile,"(a)",advance="no") "     Freq. 1     Freq. 2     Freq. 3 "
		    do idir=1,3
			    do jdir=1,3
				    do kdir=1,3
				        do ldir=1,3
                            write(ifile,"(5x,a,5x)",advance="no") dirlab(idir)//dirlab(jdir)//dirlab(kdir)//dirlab(ldir)
                        end do
                    end do
                end do
            end do
            write(ifile,*)
        end do
    end if
	
	write(*,*) "Please wait..."
	call walltime(iwalltime1)
	call fullarrange(arrg,24,4) !Generate full arrangement matrix (4!=24 :4) for gamma, each row corresponds to one permutation, e.g. 2341
! 	do i=1,24
! 		write(*,"(i4,4i1)") arrg(i,:) 
! 	end do
	do numstat=istart,iend,nstatstep
	if (isel==7) call showprog(numstat,int(dfloat(iend-1)/nstatstep)*nstatstep+1)
	do ifreq=1,nfreq !Cycle frequencies
		freq1=freqlist(ifreq,1)
		freq2=freqlist(ifreq,2)
		freq3=freqlist(ifreq,3)
		freqtot=freq1+freq2+freq3
		tmpw(1:4)=(/ -freqtot,freq1,freq2,freq3 /)
	
		do idir=1,3 !Cycle direction component
			do jdir=1,3
				do kdir=1,3
					do ldir=1,3
					
						gamma1val=0
						gamma2val=0
						tmpdir(1:4)=(/ idir,jdir,kdir,ldir /)
						do iper=1,24 !Do permutation, arrg(1,:)=1,2,3,4
							dir1=tmpdir(arrg(iper,1))
							dir2=tmpdir(arrg(iper,2))
							dir3=tmpdir(arrg(iper,3))
							dir4=tmpdir(arrg(iper,4))
							w0=tmpw(arrg(iper,1))
							w1=tmpw(arrg(iper,2))
							w2=tmpw(arrg(iper,3))
							w3=tmpw(arrg(iper,4))
							!$OMP PARALLEL SHARED(gamma1val,gamma2val) PRIVATE(istat,jstat,kstat,t1c,t2c,p1,p2,g1t,g2t) NUM_THREADS(nthreads)
							g1t=0
							g2t=0
							!$OMP DO schedule(dynamic)
							do istat=1,numstat
								do jstat=1,numstat
									!Gamma I
									do kstat=1,numstat
										t1c=trandip(istat,jstat,dir2)
										if (istat==jstat) t1c=t1c-trandip(0,0,dir2)
										t2c=trandip(jstat,kstat,dir3)
										if (jstat==kstat) t2c=t2c-trandip(0,0,dir3)
										p1=trandip(0,istat,dir1)*t1c*t2c*trandip(kstat,0,dir4)
										p2=(excene(istat)+w0)*(excene(jstat)-w2-w3)*(excene(kstat)-w3)
										g1t=g1t+p1/p2
									end do
									!Gamma II
									p1=trandip(0,istat,dir1)*trandip(istat,0,dir2)*trandip(0,jstat,dir3)*trandip(jstat,0,dir4)
									p2=(excene(istat)+w0)*(excene(istat)-w1)*(excene(jstat)-w3) !(excene(jstat)-w3) can also be (excene(jstat)+w2), they are equivalent
									g2t=g2t+p1/p2
								end do
							end do
							!$OMP END DO
							!$OMP CRITICAL
							gamma1val=gamma1val+g1t
							gamma2val=gamma2val+g2t
							!$OMP END CRITICAL
							!$OMP END PARALLEL
						end do
						gamma(idir,jdir,kdir,ldir)=gamma1val-gamma2val
						gamma1(idir,jdir,kdir,ldir)=gamma1val
                        gamma2(idir,jdir,kdir,ldir)=-gamma2val
						
					end do
				end do
			end do
		end do
		!Total Gamma
		gammaX=0;gammaY=0;gammaZ=0
		do i=1,3
			gammaX=gammaX+gamma(1,i,i,1)+gamma(1,i,1,i)+gamma(1,1,i,i)
			gammaY=gammaY+gamma(2,i,i,2)+gamma(2,i,2,i)+gamma(2,2,i,i)
			gammaZ=gammaZ+gamma(3,i,i,3)+gamma(3,i,3,i)+gamma(3,3,i,i)
		end do
		gammaX=gammaX/15;gammaY=gammaY/15;gammaZ=gammaZ/15
		gammatot=dsqrt(gammaX**2+gammaY**2+gammaZ**2)
		gammaavg1=gammaX+gammaY+gammaZ
		gammaavg2=( gamma(1,1,1,1)+gamma(2,2,2,2)+gamma(3,3,3,3) + gamma(1,1,2,2)+gamma(1,1,3,3)+gamma(2,2,3,3) &
        + gamma(2,2,1,1)+gamma(3,3,1,1)+gamma(3,3,2,2) )/5
        gamma_nor=0
        do i=1,3
            do j=1,3
                gamma_nor=gamma_nor+ 2*gamma(i,j,j,i)-gamma(i,i,j,j)
            end do
        end do
        gamma_nor=gamma_nor/15
		
		if (isel==3) then
			write(*,*) "Second hyperpolarizability tensor:"
			do jdir=1,3
				do kdir=1,3
					do ldir=1,3
						do idir=1,3
							write(*,"(2x,4a,'=',1PE14.5,2x)",advance="no") dirlab(idir),dirlab(jdir),dirlab(kdir),dirlab(ldir),gamma(idir,jdir,kdir,ldir)
							if (idir==3) write(*,*)
						end do
					end do
				end do
			end do
			write(*,"(/,' Gamma_X:',1PE14.5,'  Gamma_Y:',1PE14.5,'  Gamma_Z:',1PE14.5)") gammaX,gammaY,gammaZ
			write(*,"(a,1PE14.5)") " Magnitude of gamma:",gammatot
			write(*,"(a,1PE14.5)") " Average of gamma (definition 1), gamma ||:",gammaavg1
			write(*,"(a,1PE14.5)") " Average of gamma (definition 2):",gammaavg2
			write(*,"(a,1PE14.5)") " gamma _|_:",gamma_nor
			write(*,*)
		else if (isel==7) then
			write(10,"(i6,6(1PE14.5))") numstat,gammaX,gammaY,gammaZ,gammatot,gammaavg1,gammaavg2
            write(11,"(i10)",advance="no") numstat
            write(11,"(81(1PE14.5))",advance="no") ((((gamma(idir,jdir,kdir,ldir),ldir=1,3),kdir=1,3),jdir=1,3),idir=1,3)
            write(11,*)
            write(12,"(i10)",advance="no") numstat
            write(12,"(81(1PE14.5))",advance="no") ((((gamma1(idir,jdir,kdir,ldir),ldir=1,3),kdir=1,3),jdir=1,3),idir=1,3)
            write(12,*)
            write(13,"(i10)",advance="no") numstat
            write(13,"(81(1PE14.5))",advance="no") ((((gamma2(idir,jdir,kdir,ldir),ldir=1,3),kdir=1,3),jdir=1,3),idir=1,3)
            write(13,*)
		else if (isel==17) then
			write(10,"(3f12.6,6(1PE14.5))") freq1,freq2,freq3,gammaX,gammaY,gammaZ,gammatot,gammaavg1,gammaavg2
            write(11,"(3f12.6)",advance="no") freq1,freq2,freq3
            write(11,"(81(1PE14.5))",advance="no") ((((gamma(idir,jdir,kdir,ldir),ldir=1,3),kdir=1,3),jdir=1,3),idir=1,3)
            write(11,*)
            write(12,"(3f12.6)",advance="no") freq1,freq2,freq3
            write(12,"(81(1PE14.5))",advance="no") ((((gamma1(idir,jdir,kdir,ldir),ldir=1,3),kdir=1,3),jdir=1,3),idir=1,3)
            write(12,*)
            write(13,"(3f12.6)",advance="no") freq1,freq2,freq3
            write(13,"(81(1PE14.5))",advance="no") ((((gamma2(idir,jdir,kdir,ldir),ldir=1,3),kdir=1,3),jdir=1,3),idir=1,3)
            write(13,*)
		end if
	end do !end cycle freqlist
	end do !end cycle the number of states
	
	call walltime(iwalltime2)
	write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
	if (isel==7.or.isel==17) then
		close(10)
        close(11)
        close(12)
        close(13)
		if (isel==7) then
			write(*,*) "Done! The result has been outputted to gamma_n.txt in current folder"
			write(*,*) "The correspondence between columns and information in this file is as follows"
			write(*,*) "Column 1:  The number of states in consideration"
			write(*,*) "Column 2:  Gamma_X"
			write(*,*) "Column 3:  Gamma_Y"
			write(*,*) "Column 4:  Gamma_Z"
			write(*,*) "Column 5:  Magnitude of gamma"
			write(*,*) "Column 6:  Average of gamma (definition 1)"
			write(*,*) "Column 7:  Average of gamma (definition 2)"
            write(*,"(a)") " In addition, all components of gamma with respect to the number of considered states &
            have been exported to gamma_n_comp.txt in current folder, while its parts I and II (see gamma expression in &
            Section 3.200.8 of manual) have been respectively exported to gamma_I_n_comp.txt and gamma_II_n_comp.txt."
		else if (isel==17) then
			write(*,*) "Done! The result has been outputted to gamma_w.txt in current folder"
			write(*,*) "The correspondence between columns and information in this file is as follows"
			write(*,*) "Column 1:  Frequency of the first external field (w1)"
			write(*,*) "Column 2:  Frequency of the second external field (w2)"
			write(*,*) "Column 3:  Frequency of the third external field (w3)"
			write(*,*) "Column 4:  Gamma_X"
			write(*,*) "Column 5:  Gamma_Y"
			write(*,*) "Column 6:  Gamma_Z"
			write(*,*) "Column 7:  Magnitude of gamma"
			write(*,*) "Column 8:  Average of gamma (definition 1)"
			write(*,*) "Column 9:  Average of gamma (definition 2)"
            write(*,"(a)") " In addition, all components of gamma with respect to frequencies &
            have been exported to gamma_w_comp.txt in current folder, while its parts I and II (see gamma expression in &
            Section 3.200.8 of manual) have been respectively exported to gamma_I_w_comp.txt and gamma_II_w_comp.txt."
		end if 
	end if

! Calculate third hyperpolarizability (delta)
else if (isel==4) then
	write(*,*) "Input w1,w2,w3,w4 in a.u. for delta(-(w1+w2+w3+w4);w1,w2,w3,w4)"
	write(*,*) "e.g. 0.13,0.13,0,-0.13"
    write(*,*) "You can also input the values in nm, e.g. 450.7,-450.7,532,532 nm"
    read(*,"(a)") c80tmp
	read(c80tmp,*) freq1,freq2,freq3,freq4
	if (index(c80tmp,"nm")/=0) then
	    freq1=1239.842D0/au2eV/freq1
	    freq2=1239.842D0/au2eV/freq2
	    freq3=1239.842D0/au2eV/freq3
	    freq4=1239.842D0/au2eV/freq4
    end if
	freqtot=freq1+freq2+freq3+freq4
	tmpw(1:5)=(/ -freqtot,freq1,freq2,freq3,freq4 /)
	if (freq1/=0) then
		wavlen1=1239.842D0/(freq1*au2eV)
		write(*,"(' Wavelength of w1:',f12.6,' a.u.',f12.3,' nm')") freq1,wavlen1
    else
        write(*,*) "w1 is static"
	end if
	if (freq2/=0) then
		wavlen2=1239.842D0/(freq2*au2eV)
		write(*,"(' Wavelength of w2:',f12.6,' a.u.',f12.3,' nm')") freq2,wavlen2
    else
        write(*,*) "w2 is static"
	end if
	if (freq3/=0) then
		wavlen3=1239.842D0/(freq3*au2eV)
		write(*,"(' Wavelength of w3:',f12.6,' a.u.',f12.3,' nm')") freq3,wavlen3
    else
        write(*,*) "w3 is static"
	end if
	if (freq4/=0) then
		wavlen4=1239.842D0/(freq4*au2eV)
		write(*,"(' Wavelength of w4:',f12.6,' a.u.',f12.3,' nm')") freq4,wavlen4
    else
        write(*,*) "w4 is static"
	end if
	write(*,*)
	write(*,"(' Consider how many states? Should be <=',i6)") nstates
	read(*,*) numstat
	
	write(*,*) "Please wait patiently..."
	call walltime(iwalltime1)
	call fullarrange(arrd,120,5) !Generate full arrangement matrix (4!=24 :4) for gamma, each row corresponds to one permutation, e.g. 2341
	iprog=0
	do idir=1,3 !Cycle direction component
		do jdir=1,3
			do kdir=1,3
				do ldir=1,3
					iprog=iprog+1
					call showprog(iprog,81)
					do mdir=1,3
				
						delta1=0
						delta2=0
						delta3=0
						tmpdir(1:5)=(/ idir,jdir,kdir,ldir,mdir /)
						do iper=1,120 !Do permutation, arrd(1,:)=1,2,3,4,5
							dir1=tmpdir(arrd(iper,1))
							dir2=tmpdir(arrd(iper,2))
							dir3=tmpdir(arrd(iper,3))
							dir4=tmpdir(arrd(iper,4))
							dir5=tmpdir(arrd(iper,5))
							w0=tmpw(arrd(iper,1))
							w1=tmpw(arrd(iper,2))
							w2=tmpw(arrd(iper,3))
							w3=tmpw(arrd(iper,4))
							w4=tmpw(arrd(iper,5))
							!$OMP PARALLEL SHARED(delta1,delta2,delta3) PRIVATE(istat,jstat,kstat,lstat,t1c,t2c,t3c,p1,p2,p3,p4,d1t,d2t,d3t) NUM_THREADS(nthreads)
							d1t=0
							d2t=0
							d3t=0
							!$OMP DO schedule(dynamic)
							do istat=1,numstat
								do jstat=1,numstat
									do kstat=1,numstat
										!Delta 1
										do lstat=1,numstat
											t1c=trandip(istat,jstat,dir2)
											if (istat==jstat) t1c=t1c-trandip(0,0,dir2)
											t2c=trandip(jstat,kstat,dir3)
											if (jstat==kstat) t2c=t2c-trandip(0,0,dir3)
											t3c=trandip(kstat,lstat,dir4)
											if (kstat==lstat) t3c=t3c-trandip(0,0,dir4)
											p1=trandip(0,istat,dir1)*t1c*t2c*t3c*trandip(lstat,0,dir5)
											p2=(excene(istat)+w0)*(excene(jstat)+w0+w1)*(excene(kstat)-w3-w4)*(excene(lstat)-w4)
											d1t=d1t+p1/p2
										end do
										!Delta 2
										t3c=trandip(jstat,kstat,dir4)
										if (jstat==kstat) t3c=t3c-trandip(0,0,dir4)
										p1=trandip(0,istat,dir1)*trandip(istat,0,dir2)*trandip(0,jstat,dir3)*t3c*trandip(kstat,0,dir5)
										p2=(excene(jstat)+w2)*(excene(kstat)-w4)
										p3=1/(excene(istat)+w0)+1/(excene(istat)-w1)
										p4=1/(excene(jstat)-w3-w4)+1/(excene(kstat)+w2+w3)
										d2t=d2t+p1/p2*p3*p4
										!Delta 3, p1 is identical to delta 2 counterpart
										p2=(excene(istat)+w0)*(excene(istat)-w1)
										p3=1/(excene(jstat)-w3-w4)/(excene(kstat)-w4)
										p4=1/(excene(jstat)+w2)/(excene(kstat)+w2+w3)
										d3t=d3t+p1/p2*(p3+p4)
									end do
								end do
							end do
							!$OMP END DO
							!$OMP CRITICAL
							delta1=delta1+d1t
							delta2=delta2+d2t
							delta3=delta3+d3t
							!$OMP END CRITICAL
							!$OMP END PARALLEL
						end do
						delta(idir,jdir,kdir,ldir,mdir)=delta1-delta2/2-delta3/2
						
					end do
				end do
			end do
		end do
	end do
	
	write(*,*)
	write(*,*) "Third hyperpolarizability tensor:"
	do jdir=1,3
		do kdir=1,3
			do ldir=1,3
				do mdir=1,3
					do idir=1,3
						write(*,"(2x,5a,'=',1PE14.5,2x)",advance="no") &
						dirlab(idir),dirlab(jdir),dirlab(kdir),dirlab(ldir),dirlab(mdir),delta(idir,jdir,kdir,ldir,mdir)
						if (idir==3) write(*,*)
					end do
				end do
			end do
		end do
	end do
	call walltime(iwalltime2)
	write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
	
else if (isel==20) then !Two or three-level analysis of beta
    write(*,"(a)") " Excitation energy (E), transition electric dipole moment between ground state and excited state, &
    &and variation of dipole moment of excited states w.r.t. ground state"
    write(*,*) "                 Trans. dipole moment (a.u.)        Var. dipole moment (a.u.)"
    write(*,*) "State   E(eV)     X       Y       Z      Tot        X       Y       Z      Tot"
    do istat=1,nstates
        !write(*,"(' State',i5,'  E=',f9.4,'  X=',f10.5,' Y=',f10.5,' Z=',f10.5,' Tot=',f10.5)") &
        !istat,excene(istat)*au2eV,trandip(0,istat,:),dsqrt(sum(trandip(0,istat,:)**2))
        term1=dsqrt(sum(trandip(0,istat,1:3)**2))
        term2=dsqrt(sum(trandip(istat,istat,1:3)-trandip(0,0,1:3))**2)
        write(*,"(i5,f9.4,4f8.3,2x,4f8.3)") istat,excene(istat)*au2eV,trandip(0,istat,1:3),term1,trandip(istat,istat,1:3)-trandip(0,0,1:3),term2
    end do
    write(*,*)
    write(*,"(a)") " Input index of state(s). If only inputting one index, then two-level model analysis will be performed between ground state and this state. &
    &If inputting two indices (e.g. 2,5), three-level model analysis will be performed on them"
    write(*,"(a)") " If you want to perform two-level analysis for a range of excited states, for example 1 to 10, you should input 1-10"
    read(*,"(a)") c80tmp
    if (index(c80tmp,'-')/=0) then
        nexc=0 !Scan
        itmp=index(c80tmp,'-')
        c80tmp(itmp:itmp)=' '
        read(c80tmp,*) istat,jstat
    else if (index(c80tmp,',')/=0) then
        nexc=2
        read(c80tmp,*) istat,jstat
    else
        nexc=1
        read(c80tmp,*) istat
    end if
    if (nexc/=0) then !Two or three-level analysis for selected state(s)
        if (nexc==1) write(*,*) "Basic information of selected state:"
        if (nexc==2) write(*,*) "Basic information of selected states:"
        do i=1,nexc
            if (i==1) iexc=istat
            if (i==2) iexc=jstat
            write(*,*)
            write(*,"(' Excited state',i6)") iexc
            write(*,"(' Excitation energy',f12.6,' a.u.',f10.4,' eV')") excene(iexc),excene(iexc)*au2eV
            write(*,"(' Transition dipole moment (a.u.)')")
            write(*,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Total=',f12.6)") trandip(0,iexc,1:3),dsqrt(sum(trandip(0,iexc,1:3)**2))
            write(*,"(' Oscillator strength')") 
            write(*,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Total=',f12.6)") (2D0/3D0)*excene(iexc)*trandip(0,iexc,1:3)**2,(2D0/3D0)*excene(iexc)*sum(trandip(0,iexc,1:3)**2)
            write(*,"(' Variation of dipole moment (a.u.)')")
            write(*,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Total=',f12.6)") &
            trandip(iexc,iexc,1)-trandip(0,0,1),trandip(iexc,iexc,2)-trandip(0,0,2),trandip(iexc,iexc,3)-trandip(0,0,3),dsqrt(sum((trandip(iexc,iexc,:)-trandip(0,0,:))**2))
        end do
        if (nexc==1) then
            do idir=1,3
                term(idir)=6* (trandip(istat,istat,idir)-trandip(0,0,idir))*trandip(0,istat,idir)**2 / excene(istat)**2
            end do
            write(*,"(/,' beta evaluated by two-level model: (a.u.)')")
            write(*,"(' XXX=',f14.4,'  YYY=',f14.4,'  ZZZ=',f14.4,'  Norm=',f14.4)") term(1:3),dsqrt(sum(term(:)**2))
        else if (nexc==2) then
            term_xxx=0
            term_yyy=0
            term_zzz=0
            write(*,"(/,' Transition dipole moment between states',i6,' to',i6,': (a.u.)')") istat,jstat
            write(*,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Total=',f12.6)") trandip(istat,jstat,1:3),dsqrt(sum(trandip(istat,jstat,1:3)**2))
            term(1:3)=6* (trandip(istat,istat,1:3)-trandip(0,0,1:3))*trandip(0,istat,1:3)**2  / excene(istat)**2
            term_xxx=term_xxx+term(1)
            term_yyy=term_yyy+term(2)
            term_zzz=term_zzz+term(3)
            write(*,"(' Individual contribution of excited state',i6,' to beta: (a.u.)')") istat
            write(*,"(' XXX=',f14.4,'  YYY=',f14.4,'  ZZZ=',f14.4,'  Norm=',f14.4)") term(1:3),dsqrt(sum(term(:)**2))
            term(1:3)=6* (trandip(jstat,jstat,1:3)-trandip(0,0,1:3))*trandip(0,jstat,1:3)**2  / excene(jstat)**2
            term_xxx=term_xxx+term(1)
            term_yyy=term_yyy+term(2)
            term_zzz=term_zzz+term(3)
            write(*,"(' Individual contribution of excited state',i6,' to beta: (a.u.)')") jstat
            write(*,"(' XXX=',f14.4,'  YYY=',f14.4,'  ZZZ=',f14.4,'  Norm=',f14.4)") term(1:3),dsqrt(sum(term(:)**2))
            term(1:3)=12* trandip(istat,jstat,1:3)*trandip(0,istat,1:3)*trandip(0,jstat,1:3) / (excene(istat)*excene(jstat))
            term_xxx=term_xxx+term(1)
            term_yyy=term_yyy+term(2)
            term_zzz=term_zzz+term(3)
            write(*,"(' Coupling contribution of the two excited states to beta: (a.u.)')")
            write(*,"(' XXX=',f14.4,'  YYY=',f14.4,'  ZZZ=',f14.4,'  Norm=',f14.4)") term(1:3),dsqrt(sum(term(:)**2))
            write(*,"(/,' beta evaluated by three-level model: (a.u.)')")
            write(*,"(' XXX=',f14.4,'  YYY=',f14.4,'  ZZZ=',f14.4,'  Norm=',f14.4)") term_xxx,term_yyy,term_zzz,dsqrt(term_xxx**2+term_yyy**2+term_zzz**2)
        end if
    else !Two-level analysis for a range of excited states
        write(*,"(/,' beta evaluated by two-level model: (a.u.)')")
        do idx=istat,jstat
            do idir=1,3
                term(idir)=6* (trandip(idx,idx,idir)-trandip(0,0,idir))*trandip(0,idx,idir)**2 / excene(idx)**2
            end do
            write(*,"(' #',i5,': XXX=',f13.2,'  YYY=',f13.2,'  ZZZ=',f13.2,'  Norm=',f13.2)") idx,term(1:3),dsqrt(sum(term(:)**2))
        end do
    end if
end if

end do

end subroutine






!!---------------------------------------------------------------------------------------------
!!--------- Visualize (hyper)polarizability via unit sphere and vector representations --------
!!---------------------------------------------------------------------------------------------
subroutine vis_hypol
use util
use defvar
implicit real*8 (a-h,o-z)
real*8 :: arrowscl=0.015D0,arrowrad=0.025D0,arrowradvec=0.15D0,arrowsclvec=0.05D0
real*8 alpha(3,3),beta(3,3,3),gamma(3,3,3,3),vec(3),vectmp(3),vecrep(3),univec(3),ten2(3,3),ten3(3,3,3)
real*8,allocatable :: ptxyz(:,:),allvec(:,:)
character c200tmp*200,c80tmp*80,name*10
integer :: numpt=600,icolorarrow=1,ispecmaxlen=0

if (allocated(distmat)) then !Input file contains atom information
    sphererad=maxval(distmat)/2*b2a*1.6D0
else
    sphererad=2D0
end if

do while(.true.)
    write(*,*)
    write(*,*) "       --------------- Visualizing (hyper)polarizability ---------------"
    if (ispecmaxlen==0) write(*,*) "-8 Toggle making longest arrow on sphere has specific length, current: No"
    if (ispecmaxlen==1) write(*,*) "-8 Toggle making longest arrow on sphere has specific length, current: Yes"
    if (icolorarrow==0) write(*,*) "-7 Toggle coloring arrows on sphere, current: No"
    if (icolorarrow==1) write(*,*) "-7 Toggle coloring arrows on sphere, current: Yes"
    write(*,"(a,f9.6)") " -6 Set radius for the arrow of vector rep., current:",arrowradvec
    write(*,"(a,f9.6)") " -5 Set length scale factor for the arrow of vector rep., current:",arrowsclvec
    write(*,"(a,f9.6)") " -4 Set radius for the arrows on sphere, current:",arrowrad
    write(*,"(a,f9.6)") " -3 Set length scale factor for the arrows on sphere, current:",arrowscl
    write(*,"(a,i6)") " -2 Set number of points on the sphere, current:",numpt
    write(*,"(a,f8.3,' Angstrom')") " -1 Set sphere radius, current:",sphererad
    write(*,*) "0 Return"
    write(*,*) "1 Do analysis for polarizability (alpha)"
    write(*,*) "2 Do analysis for first-order hyperpolarizability (beta)"
    write(*,*) "3 Do analysis for second-order hyperpolarizability (gamma)"
    read(*,*) isel
    
    if (isel==-8) then
        if (ispecmaxlen==1) then
            ispecmaxlen=0
        else
            ispecmaxlen=1
        end if
    else if (isel==-7) then
        if (icolorarrow==1) then
            icolorarrow=0
        else
            icolorarrow=1
        end if
    else if (isel==-6) then
        write(*,*) "Input radius of the arrow, e.g. 0.25"
        read(*,*) arrowradvec
    else if (isel==-5) then
        write(*,*) "Input scale factor for the arrow, e.g. 0.05"
        read(*,*) arrowsclvec
    else if (isel==-4) then
        write(*,*) "Input radius of the arrows, e.g. 0.05"
        read(*,*) arrowrad
    else if (isel==-3) then
        write(*,*) "Input scale factor for the arrows, e.g. 0.15"
        read(*,*) arrowscl
    else if (isel==-2) then
        write(*,*) "Input expected number of points on the sphere, e.g. 800"
        write(*,"(a)") "  You can also input negative value to specify density of points, e.g. -8.5 means the density is 8.5 points per Angstrom^2"
        write(*,*) "Note: The actual number will be automatically maginally adjusted"
        read(*,*) numpt
        if (numpt<0) numpt=nint(4*sphererad**2*abs(numpt))
    else if (isel==-1) then
        write(*,*) "Input radius of sphere (in Angstrom), e.g. 1.5"
        read(*,*) sphererad
    else if (isel==0) then
        exit
    
    else if (isel==1.or.isel==2.or.isel==3) then
        if (isel==1) name="alpha"
        if (isel==2) name="beta"
        if (isel==3) name="gamma"
        c200tmp=trim(name)//".txt"
	    inquire(file=c200tmp,exist=alive)
        if (alive) then
            write(*,"(a)") " Since "//trim(c200tmp)//" can be found in current folder, "//trim(name)//" tensor will be directly loaded from it"
        else
            write(*,*) "Input path of the file containing full "//trim(name)//" tensor, e.g. C:\yohane.txt"
            do while(.true.)
                read(*,"(a)") c200tmp
	            inquire(file=c200tmp,exist=alive)
	            if (alive) exit
	            write(*,*) "Cannot find the file, input again!"
            end do
        end if
        write(*,*) "Loading "//trim(name)//" tensor from "//trim(c200tmp)
        open(10,file=c200tmp,status="old")
        if (isel==1) read(10,*) ((alpha(i,j),j=1,3),i=1,3)
        if (isel==2) read(10,*) (((beta(i,j,k),k=1,3),j=1,3),i=1,3)
        if (isel==3) read(10,*) ((((gamma(i,j,k,l),l=1,3),k=1,3),j=1,3),i=1,3)
        close(10)
        allocate(ptxyz(3,numpt),allvec(3,numpt))
        
        call unitspherept(ptxyz,numpt) !The inputted numpt will be automatically adjusted by this routine
        write(*,"(' Actual number of points on unit sphere:',i6)") numpt
        ptxyz=ptxyz*sphererad
        
        write(*,*) "Calculating data points on the sphere ..."
        do ipt=1,numpt
            tmpnorm=dsqrt(sum(ptxyz(:,ipt)**2))
            univec=ptxyz(:,ipt)/tmpnorm
            if (isel==1) then !alpha
                do i=1,3
                    allvec(i,ipt)=sum(alpha(i,:)*univec(:))
                end do
            else if (isel==2) then !beta
                do i=1,3
                    do j=1,3
                        ten2(i,j)=sum(beta(i,j,:)*univec(:))
                    end do
                end do
                do i=1,3
                    allvec(i,ipt)=sum(ten2(i,:)*univec(:))
                end do
            else if (isel==3) then !gamma
                do i=1,3
                    do j=1,3
                        do k=1,3
                            ten3(i,j,k)=sum(gamma(i,j,k,:)*univec(:))
                        end do
                    end do
                end do
                do i=1,3
                    do j=1,3
                        ten2(i,j)=sum(ten3(i,j,:)*univec(:))
                    end do
                end do
                do i=1,3
                    allvec(i,ipt)=sum(ten2(i,:)*univec(:))
                end do
            end if
        end do
        allvec=allvec*arrowscl
        
        valmax=0
        valmin=1E20
        do ipt=1,numpt
            tmp=dsqrt(sum(allvec(:,ipt)**2))
            if (tmp>valmax) valmax=tmp
            if (tmp<valmin) valmin=tmp
        end do
        write(*,"(' Minimal arrow length after scaling:',f12.3)") valmin
        write(*,"(' Maximal arrow length after scaling:',f12.3)") valmax
        
        if (ispecmaxlen==1) then
            write(*,*)
            write(*,*) "Input expected length of longest arrow on the sphere, e.g. 2.5"
            read(*,*) tmpmax
            sclf=tmpmax/valmax
            allvec=allvec*sclf
            valmax=valmax*sclf
            valmin=valmin*sclf
            write(*,"(' Current arrows on sphere have been further scaled by',f12.6)") sclf
        end if
        
        c80tmp=trim(name)//".tcl"
        write(*,*) "Outputting "//trim(c80tmp)//" ..."
        open(10,file=c80tmp,status="replace")
        write(10,"(a)") "color Display Background white"
        if (icolorarrow==0) then
            write(10,"(a)") "draw color white"
            do ipt=1,numpt
                call drawVMDarrow(10,ptxyz(:,ipt),allvec(:,ipt),arrowrad)
            end do
        else if (icolorarrow==1) then
            call writeVMD_BWR(10)
            do ipt=1,numpt
                tmp=dsqrt(sum(allvec(:,ipt)**2))
                !idcolor=nint(tmp/valmax*1000) !Using zero as color lower limit
                idcolor=nint((tmp-valmin)/(valmax-valmin)*1000)
                if (idcolor==0) idcolor=1
                write(10,"(a,i5)") "draw color",idcolor+50
                call drawVMDarrow(10,ptxyz(:,ipt),allvec(:,ipt),arrowrad)
            end do
        end if
        close(10)
        write(*,"(a,/)") " Done! "//trim(c80tmp)//" has been generated in current folder, it is a VMD plotting script, &
        you can run ""source "//trim(c80tmp)//""" in VMD console window to plot the map"
        deallocate(ptxyz,allvec)
        
        !Output vector representation
        vecrep=0
        if (isel==1) then
            do i=1,3
                do j=1,3
                    vecrep(i)=vecrep(i)+alpha(i,j)
                end do
            end do
            write(*,"(' Alpha_X:',1PE14.6,'   Alpha_Y:',1PE14.6,'   Alpha_Z:',1PE14.6,' a.u.')") vecrep(:)
        else if (isel==2) then
            do i=1,3
		        do j=1,3
			        vecrep(i)=vecrep(i)+beta(i,j,j)+beta(j,j,i)+beta(j,i,j)
		        end do
            end do
            vecrep=vecrep/3
            write(*,"(' Beta_X:',1PE14.6,'   Beta_Y:',1PE14.6,'   Beta_Z:',1PE14.6,' a.u.')") vecrep(:)
        else if (isel==3) then
            do i=1,3
	            do j=1,3
		            vecrep(i)=vecrep(i)+gamma(i,j,j,i)+gamma(i,j,i,j)+gamma(i,i,j,j)
	            end do
            end do
            vecrep=vecrep/15
            write(*,"(' Gamma_X:',1PE14.6,'   Gamma_Y:',1PE14.6,'   Gamma_Z:',1PE14.6,' a.u.')") vecrep(:)
        end if
        open(10,file=trim(name)//"_vec.tcl",status="replace")
        vec=0 !Origin of the arrow is (0,0,0)
        write(10,"(a)") "draw color lime"
        if (isel==2) then !beta
            call drawVMDarrow(10,vec,vecrep*arrowsclvec,arrowradvec)
        else if (isel==1.or.isel==3) then !alpha and gamma
            vectmp=0
            vectmp(1)=vecrep(1)
            call drawVMDarrow(10,vec,vectmp*arrowsclvec,arrowradvec*0.7D0)
            call drawVMDarrow(10,vec,-vectmp*arrowsclvec,arrowradvec*0.7D0)
            write(10,"(a)") "draw color pink"
            vectmp=0
            vectmp(2)=vecrep(2)
            call drawVMDarrow(10,vec,vectmp*arrowsclvec,arrowradvec*0.7D0)
            call drawVMDarrow(10,vec,-vectmp*arrowsclvec,arrowradvec*0.7D0)
            write(10,"(a)") "draw color cyan2"
            vectmp=0
            vectmp(3)=vecrep(3)
            call drawVMDarrow(10,vec,vectmp*arrowsclvec,arrowradvec*0.7D0)
            call drawVMDarrow(10,vec,-vectmp*arrowsclvec,arrowradvec*0.7D0)
        end if
        close(10)
        write(*,"(1x,a)") trim(name)//"_vec.tcl has been generated in current folder, it contains VMD command to plot "//trim(name)//" &
        tensor via vector representation"
    end if
end do

end subroutine

!!--------- Write command for drawing arrow into VMD plotting script according to inputted coordinate and arrow vector
subroutine drawVMDarrow(ifileid,ptxyz,vec,arrowrad)
integer ifileid
real*8 ptxyz(3),vec(3),arrowrad
conerad=2.5D0*arrowrad !cone radius
conepos=0.65D0
write(ifileid,"( 'draw cylinder {',3(f12.3),'} {',3(f12.3),'} radius',f5.2,' filled yes resolution 20' )") ptxyz(:),ptxyz(:)+conepos*vec(:),arrowrad
write(ifileid,"( 'draw cone {',3(f12.3),'} {',3(f12.3),'} radius',f5.2,' resolution 20' )") ptxyz(:)+conepos*vec(:),ptxyz(:)+vec(:),conerad
end subroutine

!!-------- Define 1000 customized colors (index from 51 to 1050) corresponding to variation of blue-white-red
!The reason of using 50~1050: (1) 0~32 are built-in colors (2) index >=1057 is unsupported by VMD
subroutine writeVMD_BWR(ifileid)
write(ifileid,"(a)") "set j 0                                         "
write(ifileid,"(a)") "for {set i 1} {$i<=500} {incr i} {              "
write(ifileid,"(a)") "incr j                                          "
write(ifileid,"(a)") "set red [expr double($j)/500]                       "
write(ifileid,"(a)") "set green [expr double($j)/500]                     "
write(ifileid,"(a)") "set blue 1                                      "
write(ifileid,"(a)") "color change rgb [expr $i+50] $red $green $blue"
write(ifileid,"(a)") "}                                               "
write(ifileid,"(a)") "set j 0                                         "
write(ifileid,"(a)") "for {set i 501} {$i<=1000} {incr i} {            "
write(ifileid,"(a)") "incr j                                          "
write(ifileid,"(a)") "set red 1                                       "
write(ifileid,"(a)") "set green [expr double(500-$j)/500]               "
write(ifileid,"(a)") "set blue [expr double(500-$j)/500]                "
write(ifileid,"(a)") "color change rgb [expr $i+50] $red $green $blue"
write(ifileid,"(a)") "}                                               "
end subroutine




!!-------------------------------------------------------------------
!!--------- Calculate and plot (hyper)polarizability density --------
!!-------------------------------------------------------------------
subroutine hyper_polar_dens
use defvar
use util
use functions
use GUI
implicit real*8 (a-h,o-z)
integer fieldlist(4)
character c80tmp*80,c200tmp*200,c200tmp2*200,c2000tmp*200,fname(4)*3,fpath(4)*200
real*8 :: fieldstr=0.003D0

write(*,"(/,' Note: Stepsize of finite electric field used in this module is',f7.4,' a.u., which is a good choice')") fieldstr
write(*,*)
call menutitle("Calculate and plot (hyper)polarizability density",10,1)
write(*,*) "Choose the quantity to study"
write(*,*) "0 Exit"
write(*,*) "1 Polarizability density and spatial contribution to polarizability"
write(*,*) "2 First hyperpolarizability density and spatial contribution to first hyperpolarizability"
write(*,*) "3 Second hyperpolarizability density and spatial contribution to second hyperpolarizability"
read(*,*) itype
if (itype==0) then
    return
else if (itype==1) then
    ntime=2
    fieldlist(1:2)=(/ -1,1 /)
else if (itype==2) then
    ntime=3
    fieldlist(1:3)=(/ -1,0,1 /)
else if (itype==3) then
    ntime=4
    fieldlist(1:4)=(/ -2,-1,1,2 /)
end if

write(*,*) "Choose the direction of interest"
write(*,*) "1 X"
write(*,*) "2 Y"
write(*,*) "3 Z"
read(*,*) idir

!Generate file name corresponding to different fields
if (idir==1) fname(1:ntime)(1:1)='X'
if (idir==2) fname(1:ntime)(1:1)='Y'
if (idir==3) fname(1:ntime)(1:1)='Z'
do itime=1,ntime
    if (fieldlist(itime)==0) then
        fname(itime)(2:3)="_0"
    else if (fieldlist(itime)>0) then
        write(fname(itime)(2:3),"('+',i1)") fieldlist(itime)
    else
        write(fname(itime)(2:3),"(i2)") fieldlist(itime)
    end if
end do

do while(.true.)
    write(*,"(/,a)") " 1 Generate Gaussian input files of single point task under different external electric fields"
    write(*,"(a)") " 2 Load .wfx files of single point task under different external electric fields"
    read(*,*) isel
    if (isel==1) then
        if (allocated(b)) then
            ichg=sum(a%charge)-nelec
            imult=nint(naelec)-nint(nbelec)+1
        else
            write(*,*) "Input net charge and spin multiplicity, e.g. 0,2"
            write(*,*) "If press ENTER button directly, 0,1 will be used"
            read(*,"(a)") c80tmp
            if (c80tmp==" ") then
                ichg=0
                imult=1
            else
                read(c80tmp,*) ichg,imult
            end if
        end if
        do itime=1,ntime
            c80tmp=fname(itime)//".gjf"
            write(*,*) "Generating "//trim(c80tmp)//" in current folder..."
            if (fieldlist(itime)==0) then
                c200tmp=" "
            else
                if (idir==1) c200tmp='field=X'
                if (idir==2) c200tmp='field=Y'
                if (idir==3) c200tmp='field=Z'
                !Note that definition of direction of electric field in Gaussian in inversed w.r.t. to common convention
                if (fieldlist(itime)<0) c200tmp=trim(c200tmp)//'+'
                write(c200tmp2,"(i4)") nint(-fieldstr*fieldlist(itime)*10000)
                c200tmp=trim(c200tmp)//adjustl(c200tmp2)
            end if
            open(10,file=c80tmp,status="replace")
            write(10,"(a,/)") "#P PBE1PBE/aug-cc-pVTZ int(ultrafine,acc2e=14) scf(noincfock,novaracc) out=wfx nosymm "//trim(c200tmp)
            write(10,"(a,/)") fname(itime)
            write(10,"(i2,i3)") ichg,imult
            do iatm=1,ncenter
                write(10,"(a,3f12.6)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
            end do
            write(10,*)
            write(10,"(a,/,/)") fname(itime)//".wfx"
            close(10)
        end do
        write(*,"(a)") " Generating input files finished! The computational level is set to PBE0/aug-cc-pVTZ. &
        You may manually edit the files if you want to use another level. You should run the files by Gaussian manually, &
        and then load them via option 2"
    else if (isel==2) then
        write(*,*) "Input the folder in which .wfx files calculated under different electric fields are available, e.g. D:\mari\riko\"
        write(*,*) "If press ENTER button directly, they are assumed to be in current folder"
        read(*,"(a)") c200tmp
        do itime=1,ntime
            fpath(itime)=trim(c200tmp)//fname(itime)//".wfx"
            inquire(file=fpath(itime),exist=alive)
	        if (alive) then
                write(*,"(a)") " Found "//trim(fpath(itime))
            else
	            write(*,*) "Error: Cannot find "//trim(fpath(itime))//", you need to generate it first!"
                exit
            end if
        end do
        if (itime==ntime+1) then
            write(*,*) "Well, all needed files are available"
            exit
        end if
    end if
end do

if (itype==1) then
    c200tmp="polarizability density"
    c200tmp2="spatial contribution to polarizability"
else if (itype==2) then
    c200tmp="first hyperpolarizability density"
    c200tmp2="spatial contribution to first hyperpolarizability"
else if (itype==3) then
    c200tmp="second hyperpolarizability density"
    c200tmp2="spatial contribution to second hyperpolarizability"
end if

do while(.true.)
    write(*,*)
    write(*,*) "0 Exit"
    write(*,*) "1 Generate grid data of "//trim(c200tmp)
    write(*,*) "2 Generate grid data of "//trim(c200tmp2)
    write(*,*) "3 Plot plane map of "//trim(c200tmp)
    write(*,*) "4 Plot plane map of "//trim(c200tmp2)
    read(*,*) isel
    if (isel==0) then
        return
    else if (isel==1.or.isel==2) then
        call setgrid(10,igridsel)
        if (allocated(cubmat)) deallocate(cubmat)
        if (allocated(cubmattmp)) deallocate(cubmattmp)
        allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
        !Calculate (hyper)polarizability
        do itime=1,ntime
            call dealloall(0)
            write(*,*) "Calculating "//trim(fpath(itime))
            call readinfile(fpath(itime),1)
            ifinish=0;ishowprog=1
            ntmp=floor(ny*nz/100D0)
            !$OMP PARALLEL DO SHARED(cubmattmp,ifinish,ishowprog) PRIVATE(i,j,k,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
            do k=1,nz
	            do j=1,ny
		            do i=1,nx
			            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			            cubmattmp(i,j,k)=fdens(tmpx,tmpy,tmpz)
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
            if (itype==1) then !Polarizability density
                if (itime==1) then
                    cubmat=-cubmattmp
                else if (itime==2) then
                    cubmat=cubmat+cubmattmp
                end if
            else if (itype==2) then !First hyperpolarizability density
                if (itime==1) then
                    cubmat=cubmattmp
                else if (itime==2) then
                    cubmat=cubmat-2*cubmattmp
                else if (itime==3) then
                    cubmat=cubmat+cubmattmp
                end if
            else if (itype==3) then !Second hyperpolarizability density
                if (itime==1) then
                    cubmat=-cubmattmp
                else if (itime==2) then
                    cubmat=cubmat+2*cubmattmp
                else if (itime==3) then
                    cubmat=cubmat-2*cubmattmp
                else if (itime==4) then
                    cubmat=cubmat+cubmattmp
                end if
            end if
        end do
        if (itype==1) then
            cubmat=cubmat/(2*fieldstr)
        else if (itype==2) then
            cubmat=cubmat/fieldstr**2
        else if (itype==3) then
            cubmat=cubmat/(2*fieldstr**3)
        end if
        !Transform (hyper)polarizability to spatial contribution
        if (isel==2) then
            do k=1,nz
	            do j=1,ny
		            do i=1,nx
			            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                        if (idir==1) then
                            cubmat(i,j,k)=-tmpx*cubmat(i,j,k)
                        else if (idir==2) then
                            cubmat(i,j,k)=-tmpy*cubmat(i,j,k)
                        else if (idir==3) then
                            cubmat(i,j,k)=-tmpz*cubmat(i,j,k)
                        end if
		            end do
	            end do
            end do
        end if
        deallocate(cubmattmp)
        call dealloall(0)
        write(*,*) "Reloading "//trim(firstfilename)
        call readinfile(firstfilename,1)
        call calc_dvol(dvol)
        write(*,"(/,' Integral of grid data:',f16.6,' a.u.')") sum(cubmat)*dvol
        if (isel==1) then
            sur_value=0.5D0
        else if (isel==2) then
            sur_value=2
        end if
        do while(.true.)
            write(*,*)
            write(*,*) "0 Return"
            write(*,*) "1 Visualize isosurface map"
            write(*,*) "2 Export grid data as cube file"
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==1) then
                call drawisosurgui(1)
            else if (isel2==2) then
                write(*,*) "Input path of the cube file to export, e.g. D:\Aqours\Sunshine\poldens.cub"
                write(*,*) "If press ENTER button directly, grid data will be exported to grid.cub in current folder"
                read(*,"(a)") c200tmp
                if (c200tmp==" ") c200tmp="grid.cub"
			    open(10,file=c200tmp,status="replace")
			    call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
			    close(10)
			    write(*,"(' Done! Grid data has been exported to ',a)") trim(c200tmp)
            end if
        end do
        
    else if (isel==3.or.isel==4) then
        !Generate custom operation commands
        if (allocated(custommapname)) deallocate(custommapname)
        allocate(custommapname(ntime)) !Temporarily use custommapname to pass files to be processed to study2dim
        do itime=1,ntime
            custommapname(itime)=trim(fpath(itime))
        end do
        itypetmp=itype
        if (isel==4) itypetmp=-itypetmp
        aug2D_old=aug2D
        aug2D=8
        call study2dim(2,itypetmp,idir)
        aug2D=aug2D_old
        call dealloall(0)
        write(*,*) "Reloading "//trim(firstfilename)
        call readinfile(firstfilename,1)
    end if
end do
end subroutine