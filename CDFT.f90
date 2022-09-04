!!------------ Calculate various quantities in conceptual density functional theory
subroutine CDFT
use defvar
use GUI
use functions
use util
implicit real*8 (a-h,o-z)
character keywords*200,c200tmp*200,inpname*200,selectyn,c3*3,c3p*3,c3q*3
character(len=200) wfnfile(4)
integer charge(4),spin(4) !Charge and spin multiplicity for N,N+1,N-1,N-2 states
integer :: iwcubic=0,iQCprog=1
integer :: np=1,nq=1 !Degenerate of LUMO, HOMO. Normal case is 1,1, for (quasi)-degenerate case they are p,q
real*8 atmchg(ncenter,3) !Hirshfeld charges of N,N+1,N-1,N-2 states
real*8 ene(4),E_HOMO(4) !Energy and E_HOMO of N,N+1,N-1,N-2 states
real*8 atmfneg(ncenter),atmfpos(ncenter),atmf0(ncenter),atmCDD(ncenter)
real*8 atmsneg(ncenter),atmspos(ncenter),atms0(ncenter),atmelectphi(ncenter),atmnucleophi(ncenter),atmwcubic(ncenter)
real*8 nucleophi,expterm(nmo)
real*8,allocatable :: rhoN(:,:,:),rhoNp1(:,:,:),rhoNn1(:,:,:),OW_fpos(:,:,:),OW_fneg(:,:,:)
real*8 OWfposgrid(radpot*sphpot),OWfneggrid(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),selfdens(radpot*sphpot),wfnval(nmo)
real*8,allocatable :: atmcomp(:,:),atmref(:,:)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
cubfac=1D0

write(*,"(/,a)") " !!! NOTE: If this module is used in your research, please NOT ONLY cite original paper of Multiwfn (J. Comput. Chem., 33, 580-592 (2012)), &
BUT ALSO cite the following book chapter, which comprehensively introduces feature and implementation of this module:"
write(*,"(a)") " Tian Lu, Qinxue Chen. Realization of Conceptual Density Functional Theory and Information-Theoretic Approach in &
Multiwfn Program. In Conceptual Density Functional Theory, WILEY-VCH GmbH: Weinheim (2022); pp 631-647 DOI: 10.1002/9783527829941.ch31"

do while(.true.)
    write(*,*)
    write(*,*) "---- Calculate various quantities in conceptual density functional theory ----"
    if (np==1.and.nq==1) then
		write(*,*) "-3 Set degree of degeneracy for options 1, 2 and 3, current: Nondegenerate"
    else
		write(*,"(a,i2,a,i2,a)") " -3 Set degree of degeneracy for options 1, 2 and 3, current: Degeneracy of HOMO and LUMO are",np," and",nq,", respectively"
    end if
    if (iQCprog==1) write(*,*) "-2 Choose the quantum chemistry program used in option 1, current: Gaussian"
    if (iQCprog==2) write(*,*) "-2 Choose the quantum chemistry program used in option 1, current: ORCA"
    if (np==1.and.nq==1) then
		if (iwcubic==0) write(*,*) "-1 Toggle calculating w_cubic electrophilicity index by option 2, current: No"
		if (iwcubic==1) write(*,*) "-1 Toggle calculating w_cubic electrophilicity index by option 2, current: Yes"
    end if
    write(*,*) "0 Return"
    if (iwcubic==0) write(*,"(a,i1,a,i1,a)") " 1 Generate .wfn files for N, N+",np,", N-",nq," electrons states"
    if (iwcubic==1) write(*,*) "1 Generate .wfn files for N, N+1, N-1, N-2 electrons states"
    if (iwcubic==0) write(*,*) "2 Calculate various quantitative indices"
    if (iwcubic==1) write(*,*) "2 Calculate various quantitative indices including w_cubic"
    write(*,*) "3 Calculate grid data of Fukui function and dual descriptor"
    write(*,"(a,f7.4,' a.u.')") " 4 Set delta in orbital-weighted (OW) calculation, current:",orbwei_delta
    write(*,*) "5 Print current orbital weights used in orbital-weighted (OW) calculation"
    write(*,*) "6 Calculate condensed OW Fukui function and OW dual descriptor"
    write(*,*) "7 Calculate grid data of OW Fukui function and OW dual descriptor"
    write(*,*) "8 Calculate nucleophilic and electrophilic superdelocalizabilities"
    read(*,*) isel
    
	write(c3p,"(a,i1)") "N+",np !For more convenient output
	write(c3q,"(a,i1)") "N-",nq
    
    if (isel==-3) then
        if (wfntype/=0) then
            write(*,"(a)") " Error: This function is only available for closed-shell single-determinant wavefunction! Press ENTER button to return"
            read(*,*)
            cycle
        end if
        call getHOMOidx
        idxLUMO=idxHOMO+1
        nshow=min(nmo,idxLUMO+9)-idxLUMO+1
        write(*,"(/,i3,a)") nshow," lowest unoccupied orbitals:"
        write(*,*) "E_diff denotes difference of the orbital energy with respect to LUMO"
        do imo=min(nmo,idxLUMO+9),idxLUMO,-1
            if (imo==idxLUMO) then
                write(*,"(' Orbital',i6,' (LUMO  )   Energy:',f10.3,' eV')") imo,MOene(imo)*au2eV
            else
                write(*,"(' Orbital',i6,' (LUMO+',i1,')   Energy:',f10.3,' eV  E_diff:',f10.3,' eV')") imo,imo-idxLUMO,MOene(imo)*au2eV,(MOene(imo)-MOene(idxLUMO))*au2eV
            end if
        end do
        write(*,*) "Please input degeneracy of LUMO, e.g. 3"
        read(*,*) np
        if (np>9.or.np<1) then
			if (np>9) write(*,*) "Error: The degeneracy must < 9! Press ENTER button to continue"
			if (np<1) write(*,*) "Error: The degeneracy must >=1! Press ENTER button to continue"
            read(*,*)
            np=1
            cycle
        end if
        nshow=idxHOMO-max(1,idxHOMO-9)+1
        write(*,"(/,i3,a)") nshow," highest occupied orbitals:"
        write(*,*) "E_diff denotes difference of the orbital energy with respect to HOMO"
        do imo=idxHOMO,max(1,idxHOMO-9),-1
            if (imo==idxHOMO) then
                write(*,"(' Orbital',i6,' (HOMO  )   Energy:',f10.3,' eV')") imo,MOene(imo)*au2eV
            else
                write(*,"(' Orbital',i6,' (HOMO',i2,')   Energy:',f10.3,' eV  E_diff:',f10.3,' eV')") imo,imo-idxHOMO,MOene(imo)*au2eV,(MOene(imo)-MOene(idxHOMO))*au2eV
            end if
        end do
        write(*,*) "Please input degeneracy of HOMO, e.g. 3"
        read(*,*) nq
        if (nq>9.or.nq<1) then
			if (nq>9) write(*,*) "Error: The degeneracy must < 9! Press ENTER button to continue"
			if (nq<1) write(*,*) "Error: The degeneracy must >=1! Press ENTER button to continue"
            read(*,*)
            nq=1
            cycle
        end if
		if ((np/=1.or.nq/=1).and.iwcubic==1) iwcubic=0 !w_cubic is undefined for (quasi)-degenerate case
    else if (isel==-2) then
        write(*,*) "Select the program for which the input file will be generated by option 1"
        write(*,*) "1 Gaussian"
        write(*,*) "2 ORCA"
        read(*,*) iQCprog
    else if (isel==-1) then
        if (iwcubic==1) then
            iwcubic=0
        else
            iwcubic=1
        end if
    else if (isel==2.or.isel==3) then !Check wavefunction files
        wfnfile(1)="N.wfn"
        inquire(file=wfnfile(1),exist=alive)
        if (.not.alive) then
            write(*,"(/,a)") " Unable to find N.wfn in current folder. Please input path of .wfn/wfx/fch/mwfn file of N electrons state, e.g. /sob/N.fch (Note that molden file should not be used, as it does not record system energy)"
            read(*,"(a)") c200tmp
            inquire(file=c200tmp,exist=alive)
            if (.not.alive) then
                write(*,*) "Error: Unable to find this file!"
                cycle
            end if
            wfnfile(1)=c200tmp
        else
            write(*,*) "N.wfn has been found in current folder"
        end if
        wfnfile(2)=c3p//".wfn"
        inquire(file=wfnfile(2),exist=alive)
        if (.not.alive) then
            write(*,"(/,a)") " Unable to find "//c3p//".wfn in current folder. Please input path of .wfn/wfx/fch/mwfn file of "//c3p//" electrons state, e.g. /sob/"//c3p//".fch"
            read(*,"(a)") c200tmp
            inquire(file=c200tmp,exist=alive)
            if (.not.alive) then
                write(*,*) "Error: Unable to find this file!"
                cycle
            end if
            wfnfile(2)=c200tmp
        else
            write(*,*) c3p//".wfn has been found in current folder"
        end if
        wfnfile(3)=c3q//".wfn"
        inquire(file=wfnfile(3),exist=alive)
        if (.not.alive) then
            write(*,"(/,a)") " Unable to find "//c3q//".wfn in current folder. Please input path of .wfn/wfx/fch/mwfn file of "//c3q//" electrons state, e.g. /sob/"//c3q//".fch"
            read(*,"(a)") c200tmp
            inquire(file=c200tmp,exist=alive)
            if (.not.alive) then
                write(*,*) "Error: Unable to find this file!"
                cycle
            end if
            wfnfile(3)=c200tmp
        else
            write(*,*) c3q//".wfn has been found in current folder"
        end if
        if (iwcubic==1) then
            wfnfile(4)="N-2.wfn"
            inquire(file=wfnfile(4),exist=alive)
            if (.not.alive) then
                write(*,"(/,a)") " Unable to find N-2.wfn in current folder. Please input path of .wfn/wfx/fch/mwfn file of N-2 electrons state, e.g. /sob/N-2.fch"
                read(*,"(a)") c200tmp
                inquire(file=c200tmp,exist=alive)
                if (.not.alive) then
                    write(*,*) "Error: Unable to find this file!"
                    cycle
                end if
                wfnfile(4)=c200tmp
            else
                write(*,*) "N-2.wfn has been found in current folder"
            end if
        end if
    end if
    
    if (isel==0) then
        return
        
    else if (isel==1) then
        if (iQCprog==1) then
            write(*,*) "Input Gaussian keywords used for single point task, e.g. PBE1PBE/def2SVP"
            write(*,*) "You can also meantime add some keywords for facilitating SCF convergence"
            write(*,*) "If press ENTER button directly, B3LYP/6-31G* will be employed"
            read(*,"(a)") keywords
            if (keywords==" ") keywords="B3LYP/6-31G* SCF=conver=6"
            c200tmp=keywords
            if (index(keywords,"gen")/=0) then
                inquire(file="basis.txt",exist=alive)
                if (.not.alive) then
                    write(*,*) "Error: ""gen"" keyword is found, but basis.txt cannot be found in current folder!"
                    write(*,*) "Press ENTER button to exit"
                    read(*,*)
                    exit
                end if
            end if
            !"nosymm" is not absolutely needed, however, because Gaussian may move the coordinate, the final coordinate in the resulting .wfn
            !may be different to the coordinate in the firstly loaded file, therefore add nosymm to guarantee the coordinate consistency
            keywords="#P "//trim(c200tmp)//" out=wfn nosymm"
        
            write(*,*) "Input the net charge and spin multiplicity for N electrons state, e.g. 0 1"
            if (np==1.and.nq==1) then !Normal case
				if (iwcubic==0) write(*,"(a)") " Note: If pressing ENTER button directly, (0 1), (-1 2) and (1 2) will be employed for N, N+1 and N-1 states, respectively"
				if (iwcubic==1) write(*,"(a)") " Note: If pressing ENTER button directly, (0 1), (-1 2), (1 2) and (2,1) will be employed for N, N+1, N-1 and N-2 states, respectively"
            else
				write(*,"(a,i1,1x,i1,a,i1,1x,i1,a)") " Note: If pressing ENTER button directly, (0 1), (-",np,np+1,") and (",nq,nq+1,") will be employed for N, "//c3p//" and "//c3q//" states, respectively"
            end if
            read(*,"(a)") c200tmp
            if (c200tmp==" ") then
				if (np==1.and.nq==1) then !Normal case
					charge(1)=0;spin(1)=1
					charge(2)=-1;spin(2)=2
					charge(3)=1;spin(3)=2
					charge(4)=2;spin(4)=1
                else
					charge(1)=0;spin(1)=1
					charge(2)=-np;spin(2)=np+1
					charge(3)=nq;spin(3)=nq+1
                end if
            else
                read(c200tmp,*) charge(1),spin(1)
                write(*,"(a,i1,a,2i2)") " Input the net charge and spin multiplicity for N+",np," electrons state, e.g. ",-np,np+1
                read(*,*) charge(2),spin(2)
                write(*,"(a,i1,a,2i2)") " Input the net charge and spin multiplicity for N-",nq," electrons state, e.g. ",nq,nq+1
                read(*,*) charge(3),spin(3)
                if (iwcubic==1) then
                    write(*,*) "Input the net charge and spin multiplicity for N-2 electrons state, e.g. 2 1"
                    read(*,*) charge(4),spin(4)
                end if
            end if
        
            !Generate .gjf for N, N+1, N-1, N-2 states
            nstates=3
            if (iwcubic==1) nstates=4
            do istate=1,nstates
                if (istate==1) inpname="N.gjf"
                if (istate==2) inpname=c3p//".gjf"
                if (istate==3) inpname=c3q//".gjf"
                if (istate==4) inpname="N-2.gjf"
                write(*,*) "Generating "//trim(inpname)//"..."
        
                open(10,file=inpname)
                !Under windows, most user don't know Default.Rou must be placed in current folder to make it take effect. So directly set number of cores
                if (isys==1)  write(10,"(a,i4)") "%nprocs=",nthreads
                write(10,"(a)") trim(keywords)
                write(10,*)
                write(10,"(a)") "Generated by Multiwfn"
                write(10,*)
                write(10,"(2i2)") charge(istate),spin(istate)
                do iatm=1,ncenter
	                write(10,"(a,1x,3f14.8)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
                end do
                if (index(keywords,"gen")/=0) then
                    open(11,file="basis.txt")
                    write(10,*)
                    do while(.true.)
                        read(11,"(a)",iostat=ierror) c200tmp
                        if (ierror/=0.or.c200tmp==" ") exit
                        write(10,"(a)") trim(c200tmp)
                    end do
                    close(11)
                end if
                write(10,*)
                if (istate==1) write(10,"(a)") "N.wfn"
                if (istate==2) write(10,"(a)") c3p//".wfn"
                if (istate==3) write(10,"(a)") c3q//".wfn"
                if (istate==4) write(10,"(a)") "N-2.wfn"
                write(10,*)
                write(10,*)
                close(10)
            end do
            write(*,*) "Gaussian input files for all states have been generated in current folder"
            write(*,*)
            inquire(file=gaupath,exist=alive)
            selectyn='n'
            if (alive) then
                write(*,"(a)") " Do you want to invoke Gaussian to calculate these .gjf files now to yield .wfn &
                files, and then automatically delete the .gjf and .out files? (y/n)"
                write(*,*) "Note: You can manually edit the .gjf files before inputting ""y"""
                read(*,*) selectyn
                if (selectyn=='y'.or.selectyn=='Y') then
                    call runGaussian("N.gjf",isuccess1)
                    if (isuccess1==1) then
                        call delfile("N.gjf N.out")
                        write(*,"(a,/)") " Now current folder should contain N.wfn"
                    else
                        write(*,"(a,/)") " The task has failed! Please manually check N.gjf and N.out"
                    end if
                    call runGaussian(c3p//".gjf",isuccess2)
                    if (isuccess2==1) then
                        call delfile(c3p//".gjf "//c3p//".out")
                        write(*,"(a,/)") " Now current folder should contain "//c3p//".wfn"
                    else
                        write(*,"(a,/)") " The task has failed! Please manually check "//c3p//".gjf and "//c3p//".out"
                    end if
                    call runGaussian(c3q//".gjf",isuccess3)
                    if (isuccess3==1) then
                        call delfile(c3q//".gjf "//c3q//".out")
                        write(*,"(a,/)") " Now current folder should contain "//c3q//".wfn"
                    else
                        write(*,"(a,/)") " The task has failed! Please manually check "//c3q//".gjf and "//c3q//".out"
                    end if
                    if (iwcubic==1) then
                        call runGaussian("N-2.gjf",isuccess4)
                        if (isuccess4==1) then
                            call delfile("N-2.gjf N-2.out")
                            write(*,"(a,/)") " Now current folder should contain N-2.wfn"
                        else
                            write(*,"(a,/)") " The task has failed! Please manually check N-2.gjf and N-2.out"
                        end if
                    end if
                    if (iwcubic==0.and.isuccess1*isuccess2*isuccess3==1) then
                        write(*,"(a)") " Since N.wfn, "//c3p//".wfn and "//c3q//".wfn have been successfully generated, &
                        now you can use option 2 or 3 to start the analysis"
                    else if (iwcubic==1.and.isuccess1*isuccess2*isuccess3*isuccess4==1) then
                        write(*,"(a)") " Since N.wfn, N+1.wfn, N-1.wfn and N-2.wfn have been successfully generated, &
                        now you can use option 2 or 3 to start the analysis"
                    else
                        write(*,"(a)") " Since one or more .wfn file was not successfully generated, the analysis cannot be conducted currently"
                    end if
                end if
            else
                write(*,"(a)") " Since ""gaupath"" in settings.ini has not been set to actual path of Gaussian executable, &
                automatically invoking Gaussian to run the input files is skipped"
                write(*,*)
            end if
            if ((.not.alive).or.selectyn=='n') then
				write(*,"(a)") " Now please manually run the input files by Gaussian, and then put the generated .wfn files to current folder, so that &
				options 2 and 3 can perform analyses based on them"
            end if
        else if (iQCprog==2) then !ORCA
            write(*,*) "Input ORCA keywords used single point task, e.g. PBE0 def2-TZVP"
            write(*,*) "You can also meantime add some keywords for facilitating SCF convergence"
            write(*,*) "If press ENTER button directly, B3LYP/G 6-31G* will be employed"
            read(*,"(a)") keywords
            if (keywords==" ") keywords="B3LYP/G 6-31G* autoaux"
        
            write(*,*) "Input the net charge and spin multiplicity for N electrons state, e.g. 0 1"
            if (np==1.and.nq==1) then
				if (iwcubic==0) write(*,"(a)") " Note: If pressing ENTER button directly, (0 1), (-1 2) and (1 2) will be employed for N, N+1 and N-1 states, respectively"
				if (iwcubic==1) write(*,"(a)") " Note: If pressing ENTER button directly, (0 1), (-1 2), (1 2) and (2,1) will be employed for N, N+1, N-1 and N-2 states, respectively"
            else
				write(*,"(a,i1,1x,i1,a,i1,1x,i1,a)") " Note: If pressing ENTER button directly, (0 1), (-",np,np+1,") and (",nq,nq+1,") will be employed for N, "//c3p//" and "//c3q//" states, respectively"
            end if
            read(*,"(a)") c200tmp
            if (c200tmp==" ") then
				if (np==1.and.nq==1) then !Normal case
					charge(1)=0;spin(1)=1
					charge(2)=-1;spin(2)=2
					charge(3)=1;spin(3)=2
					charge(4)=2;spin(4)=1
                else
					charge(1)=0;spin(1)=1
					charge(2)=-np;spin(2)=np+1
					charge(3)=nq;spin(3)=nq+1
                end if
            else
                read(c200tmp,*) charge(1),spin(1)
                write(*,"(a,i1,a,2i2)") " Input the net charge and spin multiplicity for N+",np," electrons state, e.g. ",-np,np+1
                read(*,*) charge(2),spin(2)
                write(*,"(a,i1,a,2i2)") " Input the net charge and spin multiplicity for N-",nq," electrons state, e.g. ",nq,nq+1
                read(*,*) charge(3),spin(3)
                if (iwcubic==1) then
                    write(*,*) "Input the net charge and spin multiplicity for N-2 electrons state, e.g. 2 1"
                    read(*,*) charge(4),spin(4)
                end if
            end if
        
            !Generate .inp for N, N+1, N-1, N-2 states
            nstates=3
            if (iwcubic==1) nstates=4
            do istate=1,nstates
                if (istate==1) inpname="N.inp"
                if (istate==2) inpname=c3p//".inp"
                if (istate==3) inpname=c3q//".inp"
                if (istate==4) inpname="N-2.inp"
                write(*,*) "Generating "//trim(inpname)//"..."
                open(10,file=inpname)
                write(10,"(a)") "! "//trim(keywords)//" aim"
                write(10,"('%pal nprocs',i4,' end')") nthreads
                write(10,"(a,i6)") "%maxcore",1000
                write(10,"('* xyz ',2i3)") charge(istate),spin(istate)
                do iatm=1,ncenter
	                write(10,"(a,1x,3f14.8)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
                end do
                write(10,"(a)") "*"
                close(10)
            end do
            
            write(*,*) "ORCA input files for all states have been generated in current folder"
            inquire(file=orcapath,exist=alive)
            selectyn='n'
            if (alive) then
                write(*,*)
                write(*,"(a)") " Do you want to invoke ORCA to calculate these .inp files now to yield .wfn &
                files, and then automatically delete the .inp, .out and temporary files? (y/n)"
                write(*,*) "Note: You can manually edit the .inp files before inputting ""y"""
                read(*,*) selectyn
                if (selectyn=='y'.or.selectyn=='Y') then
                    call runORCA("N.inp",isuccess1)
                    if (isuccess1==1) then
                        call delfile("N.inp N.out")
                        write(*,"(a)") " Now current folder should contain N.wfn"
                    else
                        write(*,"(a)") " The task has failed! Please manually check N.inp and N.out"
                    end if
                    call delfile("N.gbw N.densities N.PDAT.tmp N_property.txt N.wfx")
                    write(*,*)
                    call runORCA(c3p//".inp",isuccess2)
                    if (isuccess2==1) then
                        call delfile(c3p//".inp "//c3p//".out")
                        write(*,"(a)") " Now current folder should contain "//c3p//".wfn"
                    else
                        write(*,"(a)") " The task has failed! Please manually check "//c3p//".inp and "//c3p//".out"
                    end if
                    call delfile(c3p//".gbw "//c3p//".densities "//c3p//".PDAT.tmp "//c3p//"_property.txt "//c3p//".wfx")
                    write(*,*)
                    call runORCA(c3q//".inp",isuccess3)
                    if (isuccess3==1) then
                        call delfile(c3q//".inp "//c3q//".out")
                        write(*,"(a)") " Now current folder should contain "//c3q//".wfn"
                    else
                        write(*,"(a)") " The task has failed! Please manually check "//c3q//".inp and "//c3q//".out"
                    end if
                    call delfile(c3q//".gbw "//c3q//".densities "//c3q//".PDAT.tmp "//c3q//"_property.txt "//c3q//".wfx")
                    write(*,*)
                    if (iwcubic==1) then
                        call runORCA("N-2.inp",isuccess4)
                        if (isuccess4==1) then
                            call delfile("N-2.inp N-2.out")
                            write(*,"(a)") " Now current folder should contain N-2.wfn"
                        else
                            write(*,"(a)") " The task has failed! Please manually check N-2.inp and N-2.out"
                        end if
                        call delfile("N-2.gbw N-2.densities N-2.PDAT.tmp N-2_property.txt N-2.wfx")
                        write(*,*)
                    end if
                    if (iwcubic==0.and.isuccess1*isuccess2*isuccess3==1) then
                        write(*,"(a)") " Since N.wfn, "//c3p//".wfn and "//c3q//".wfn have been successfully generated, &
                        now you can use option 2 or 3 to start the analysis"
                    else if (iwcubic==1.and.isuccess1*isuccess2*isuccess3*isuccess4==1) then
                        write(*,"(a)") " Since N.wfn, N+1.wfn, N-1.wfn and N-2.wfn have been successfully generated, &
                        now you can use option 2 or 3 to start the analysis"
                    else
                        write(*,"(a)") " Since one or more .wfn file was not successfully generated, the analysis cannot be conducted currently"
                    end if
                end if
            else
                write(*,"(a)") " Since ""orcapath"" in settings.ini has not been set to actual path of ORCA executable, &
                automatically invoking ORCA to run the input files is skipped"
                write(*,*)
            end if
            if ((.not.alive).or.selectyn=='n') then
				write(*,"(a)") " Now please manually run the input files by ORCA, and then put the generated .wfn files to current folder, so that &
				options 2 and 3 can perform analyses based on them"
            end if
        end if
    
    else if (isel==2) then
        write(*,"(/,' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot
        call dealloall(0)
        !N electrons
        call readinfile(wfnfile(1),1)
        ene(1)=totenergy
        call getHOMOidx
        E_HOMO(1)=max(MOene(idxHOMO),MOene(idxHOMOb))
        write(*,*) "Calculating Hirshfeld charges for N electrons state..."
        call genHirshfeld(atmchg(:,1))
        call dealloall(0)
        !N+1 electrons
        call readinfile(wfnfile(2),1)
        ene(2)=totenergy
        call getHOMOidx
        E_HOMO(2)=max(MOene(idxHOMO),MOene(idxHOMOb))
        write(*,"(a,i1,a)") " Calculating Hirshfeld charges for N+",np," electrons state..."
        call genHirshfeld(atmchg(:,2))
        call dealloall(0)
        !N-1 electrons
        call readinfile(wfnfile(3),1)
        ene(3)=totenergy
        call getHOMOidx
        E_HOMO(3)=max(MOene(idxHOMO),MOene(idxHOMOb))
        write(*,"(a,i1,a)") " Calculating Hirshfeld charges for N-",nq," electrons state..."
        call genHirshfeld(atmchg(:,3))
        call dealloall(0)
        !N-2 electrons
        if (iwcubic==1) then
            call readinfile(wfnfile(4),1)
            ene(4)=totenergy
            call getHOMOidx
            E_HOMO(4)=max(MOene(idxHOMO),MOene(idxHOMOb))
            call dealloall(0)
        end if
        write(*,*) "Reloading the file initially loaded after booting up Multiwfn..."
        call readinfile(firstfilename,1)
        write(*,*)
        
        VIP=(ene(3)-ene(1))/nq ![E(N-q) - E(N)]/nq
        VEA=(ene(1)-ene(2))/np ![E(N) - E(N+p)]/np
        elenegMul=(VIP+VEA)/2D0
        chempot=-elenegMul
        hardness=VIP-VEA
        softness=1/hardness
        electphi=chempot**2/(2*hardness)
        E_HOMO_TCE=-0.335198D0 !TCE at B3LYP/6-31G* wavefunction and geometry
        nucleophi=E_HOMO(1)-E_HOMO_TCE
        
        atmfneg(:)=(atmchg(:,3)-atmchg(:,1))/nq
        atmfpos(:)=(atmchg(:,1)-atmchg(:,2))/np
        atmCDD(:)=atmfpos(:)-atmfneg(:)
        atmf0(:)=(atmfpos(:)+atmfneg(:))/2
        atmsneg(:)=atmfneg(:)*softness
        atmspos(:)=atmfpos(:)*softness
        atms0(:)=atmf0(:)*softness
        atmelectphi(:)=atmfpos(:)*electphi !e*Hartree
        atmnucleophi(:)=atmfneg(:)*nucleophi !e*Hartree
        
        if (iwcubic==1) then
            !Calculate terms specific for w_cubic
            VIP2=ene(4)-ene(3) !E(N-2) - E(N-1)
            c_miu=(-2*VEA-5*VIP+VIP2)/6
            c_eta=hardness
            c_gamma=2*VIP-VIP2-VEA
            w_cubic=c_miu**2/(2*c_eta)*(1+c_miu*c_gamma/3/c_eta**2)
            atmwcubic(:)=atmfpos(:)*w_cubic !e*Hartree
        end if
        
        open(10,file="CDFT.txt",status="replace")
        if (wfntype==0.or.wfntype==1.or.wfntype==2) then
            write(10,"(a)") " Note: the E(HOMO) of TCE used for evaluating nucleophilicity index is the value evaluated at B3LYP/6-31G* level"
        else
            write(10,"(a)") " Note: Nucleophilicity index was not calculated because E(HOMO) is not available at present theoretical method"
        end if
        write(10,*)
        write(10,*) "Hirshfeld charges, condensed Fukui functions and condensed dual descriptors"
        write(10,*) "Units used below are ""e"" (elementary charge)"
        write(10,*) "    Atom     q(N)    q("//c3p//")   q("//c3q//")     f-       f+       f0      CDD"
        do iatm=1,ncenter
            write(10,"(i6,'(',a')',7f9.4)") iatm,a(iatm)%name,atmchg(iatm,:),atmfneg(iatm),atmfpos(iatm),atmf0(iatm),atmCDD(iatm)
        end do
        write(10,*)
        write(10,*) "Condensed local electrophilicity/nucleophilicity index (e*eV)"
        write(10,*) "    Atom              Electrophilicity          Nucleophilicity"
        do iatm=1,ncenter
            write(10,"(i6,'(',a')',2f25.5)") iatm,a(iatm)%name,atmelectphi(iatm)*au2eV,atmnucleophi(iatm)*au2eV
        end do
        if (iwcubic==1) then
            write(10,*)
            write(10,*) "Condensed local cubic electrophilicity index (e*eV)"
            write(10,*) "    Atom              Value"
            do iatm=1,ncenter
                write(10,"(i6,'(',a')',f19.5)") iatm,a(iatm)%name,atmwcubic(iatm)*au2eV
            end do
        end if
        write(10,*)
        write(10,"(a)") " Condensed local softnesses (Hartree*e) and relative electrophilicity/nucleophilicity (dimensionless)"
        write(10,*) "    Atom         s-          s+          s0        s+/s-       s-/s+"
        do iatm=1,ncenter
            write(10,"(i6,'(',a')',5f12.4)") iatm,a(iatm)%name,atmsneg(iatm),&
            atmspos(iatm),atms0(iatm),atmspos(iatm)/atmsneg(iatm),atmsneg(iatm)/atmspos(iatm)
        end do
        
        write(10,*)
        write(10,"(a,f14.6,' Hartree')") " E(N):  ",ene(1)
        write(10,"(a,f14.6,' Hartree')") " E("//c3p//"):",ene(2)
        write(10,"(a,f14.6,' Hartree')") " E("//c3q//"):",ene(3)
        if (wfntype==0.or.wfntype==1.or.wfntype==2) then
            write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " E_HOMO(N):  ",E_HOMO(1),E_HOMO(1)*au2eV
            write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " E_HOMO("//c3p//"):",E_HOMO(2),E_HOMO(2)*au2eV
            write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " E_HOMO("//c3q//"):",E_HOMO(3),E_HOMO(3)*au2eV
        end if
        write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Vertical IP:",VIP,VIP*au2eV
        if (iwcubic==1) write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Vertical second IP:",VIP2,VIP2*au2eV
        write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Vertical EA:",VEA,VEA*au2eV
        write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Mulliken electronegativity: ",elenegMul,elenegMul*au2eV
        write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Chemical potential:         ",chempot,chempot*au2eV
        write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Hardness (=fundamental gap):",hardness,hardness*au2eV
        write(10,"(a,f12.6,' Hartree^-1,',f10.4,' eV^-1')") " Softness:",softness,softness/au2eV
        write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Electrophilicity index:",electphi,electphi*au2eV
       if (wfntype==0.or.wfntype==1.or.wfntype==2) write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Nucleophilicity index: ",nucleophi,nucleophi*au2eV
        if (iwcubic==1) write(10,"(a,f12.6,' Hartree,',f10.4,' eV')") " Cubic electrophilicity index (w_cubic):",w_cubic,w_cubic*au2eV
        close(10)
        write(*,*) "Done! Data have been outputted to CDFT.txt in current folder!"
        
    else if (isel==3) then
        if (allocated(cubmat)) deallocate(cubmat)
        aug3Dold=aug3D
        aug3D=3 !Commonly 3 Bohr extension distance is adequate
        write(*,*)
        call setgrid(1,inouse)
        aug3D=aug3Dold
        allocate(rhoN(nx,ny,nz),rhoNp1(nx,ny,nz),rhoNn1(nx,ny,nz),cubmat(nx,ny,nz))
        call dealloall(0)
        !N electrons
        call readinfile(wfnfile(1),1)
        allocate(atmref(3,ncenter)) !Used to check consistency of coordinate between different wavefunction files
        atmref(1,:)=a%x;atmref(2,:)=a%y;atmref(3,:)=a%z
        write(*,*) "Calculating electron density grid data for N electrons state..."
        call savecubmat(1,1,0)
        rhoN=cubmat
        call dealloall(0)
        !N+1 electrons
        call readinfile(wfnfile(2),1)
        do iatm=1,ncenter
            devdist=dsqrt((atmref(1,iatm)-a(iatm)%x)**2+(atmref(2,iatm)-a(iatm)%y)**2+(atmref(3,iatm)-a(iatm)%z)**2)*b2a
            if (devdist>0.05D0) then
                write(*,"(' Warning: Coodinate of atom',i5,' in N+1 state file is different to that in N state file!')") iatm
                write(*,"(' X,Y,Z of the atom in N state file:  ',3f11.5,' Angstrom')") atmref(:,iatm)*b2a
                write(*,"(' X,Y,Z of the atom in N+1 state file:',3f11.5,' Angstrom')") a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%x*b2a
                write(*,"(a)") " You need to carefully check input files of quantum chemistry program for generating the wavefunction files so that all &
                atomic coordinates in these files are completely identical, otherwise the resulting density difference will be fully meaningless!"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                exit
            end if
        end do
        write(*,*) "Calculating electron density grid data for "//c3p//" electrons state..."
        call savecubmat(1,1,0)
        rhoNp1=cubmat
        call dealloall(0)
        !N-1 electrons
        call readinfile(wfnfile(3),1)
        do iatm=1,ncenter
            devdist=dsqrt((atmref(1,iatm)-a(iatm)%x)**2+(atmref(2,iatm)-a(iatm)%y)**2+(atmref(3,iatm)-a(iatm)%z)**2)*b2a
            if (devdist>0.05D0) then
                write(*,"(' Warning: Coodinate of atom',i5,' in N-1 state file is different to that in N state file!')") iatm
                write(*,"(' X,Y,Z of the atom in N state file:  ',3f11.5,' Angstrom')") atmref(:,iatm)*b2a
                write(*,"(' X,Y,Z of the atom in N-1 state file:',3f11.5,' Angstrom')") a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%x*b2a
                write(*,"(a)") " You need to carefully check input files of quantum chemistry program for generating the wavefunction files so that all &
                atomic coordinates in these files are completely identical, otherwise the resulting density difference will be fully meaningless!"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                exit
            end if
        end do
        write(*,*) "Calculating electron density grid data for "//c3q//" electrons state..."
        call savecubmat(1,1,0)
        rhoNn1=cubmat
        call dealloall(0)
        deallocate(atmref)
        write(*,*) "Loading the file initially loaded after booting up Multiwfn..."
        call readinfile(firstfilename,1)
        
	 	sur_value=0.01D0
        do while(.true.)
            write(*,*)
            write(*,"(' -1 Set the value multiplied to all grid data, current:',f12.6)") cubfac
            write(*,*) "0 Return"
            write(*,*) "1 Visualize isosurface of f+"
            write(*,*) "2 Visualize isosurface of f-"
            write(*,*) "3 Visualize isosurface of f0"
            write(*,*) "4 Visualize isosurface of dual descriptor"
            write(*,*) "5 Export grid data of f+ as f+.cub in current folder"
            write(*,*) "6 Export grid data of f- as f-.cub in current folder"
            write(*,*) "7 Export grid data of f0 as f0.cub in current folder"
            write(*,*) "8 Export grid data of dual descriptor as DD.cub in current folder"
            read(*,*) isel2
            if (isel2==-1) then
                write(*,*) "Input the value, e.g. 0.325"
                read(*,*) cubfac
            else if (isel2==0) then
                deallocate(rhoN,rhoNp1,rhoNn1)
                exit
            else
                if (isel2==1.or.isel2==5) cubmat=(rhoNp1-rhoN)/np !f+
                if (isel2==2.or.isel2==6) cubmat=(rhoN-rhoNn1)/nq !f-
                if (isel2==3.or.isel2==7) cubmat=((rhoNp1-rhoN)/np+(rhoN-rhoNn1)/nq)/2 !(f+ + f-)/2
                if (isel2==4.or.isel2==8) cubmat=(rhoNp1-rhoN)/np - (rhoN-rhoNn1)/nq !f+ - f-
                cubmat=cubmat*cubfac
                if (isel2<=4) then
		            call drawisosurgui(1)
                else
                    if (isel2==5) open(10,file="f+.cub",status="replace")
                    if (isel2==6) open(10,file="f-.cub",status="replace")
                    if (isel2==7) open(10,file="f0.cub",status="replace")
                    if (isel2==8) open(10,file="DD.cub",status="replace")
                    call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
                    close(10)
                    write(*,*) "Exporting finished!"
                end if
            end if
        end do
    end if
    
    if (isel==5.or.isel==6.or.isel==7.or.isel==8) then
        if (.not.allocated(CObasa)) then
            write(*,"(a)") " Error: Only mwfn/fch/molden/gms file can be used for this function. Press ENTER button to return"
            read(*,*)
            cycle
        end if
        if (wfntype/=0) then
            write(*,"(a)") " Error: Only closed-shell single-determinant wavefunction is supported by this function. Press ENTER button to return"
            read(*,*)
            cycle
        end if
    end if
    
    !Calculate orbital-weighted functions
    if (isel==4) then
        write(*,*) "Input delta parameter, e.g. 0.015 (in Hartree)"
        read(*,*) orbwei_delta
        
    else if (isel==5.or.isel==6) then
        call getHOMOidx
        idxLUMO=idxHOMO+1
        chempot=(MOene(idxHOMO)+MOene(idxLUMO))/2
        do imo=1,nmo
            expterm(imo)=exp( -((chempot-MOene(imo))/orbwei_delta)**2 )
        end do
        denomin_pos=sum(expterm(idxLUMO:nmo))
        denomin_neg=sum(expterm(1:idxHOMO))
        if (isel==5) then !Show orbital weights
            write(*,"(' HOMO energy:',f12.6,' a.u.',f12.3,' eV')") MOene(idxHOMO),MOene(idxHOMO)*au2eV
            write(*,"(' LUMO energy:',f12.6,' a.u.',f12.3,' eV')") MOene(idxLUMO),MOene(idxLUMO)*au2eV
            write(*,"(' Chemical potential:',f12.6,' a.u.',f12.3,' eV')") chempot,chempot*au2eV
            write(*,"(' Delta parameter:   ',f12.6,' a.u.',f12.3,' eV')") orbwei_delta,orbwei_delta*au2eV
            write(*,"(a)") " The ""E_diff"" in below output denotes difference between orbital energy with respect to chemical potential"
            write(*,*)
            write(*,*) "10 Highest weights in orbital-weighted f+"
            totwei=0
            do imo=idxLUMO,min(nmo,idxHOMO+9)
                weithis=expterm(imo)/denomin_pos
                if (imo==idxLUMO) then
                    write(*,"(' Orbital',i6,' (LUMO  )   Weight:',f7.2,' %   E_diff:',f10.3,' eV')") imo,weithis*100,(MOene(imo)-chempot)*au2eV
                else
                    write(*,"(' Orbital',i6,' (LUMO+',i1,')   Weight:',f7.2,' %   E_diff:',f10.3,' eV')") imo,imo-idxLUMO,weithis*100,(MOene(imo)-chempot)*au2eV
                end if
                totwei=totwei+weithis
            end do
            write(*,"(' Total weight of above listed orbitals:',f7.2,' %')") totwei*100
            write(*,*)
            write(*,*) "10 Highest weights in orbital-weighted f-"
            totwei=0
            do imo=idxHOMO,max(1,idxHOMO-9),-1
                weithis=expterm(imo)/denomin_neg
                if (imo==idxHOMO) then
                    write(*,"(' Orbital',i6,' (HOMO  )   Weight:',f7.2,' %   E_diff:',f10.3,' eV')") imo,weithis*100,(MOene(imo)-chempot)*au2eV
                else
                    write(*,"(' Orbital',i6,' (HOMO',i2,')   Weight:',f7.2,' %   E_diff:',f10.3,' eV')") imo,imo-idxHOMO,weithis*100,(MOene(imo)-chempot)*au2eV
                end if
                totwei=totwei+weithis
            end do
            write(*,"(' Total weight of above listed orbitals:',f7.2,' %')") totwei*100
        else if (isel==6) then !Calculate condensed values
            if (iautointgrid==1) then
                ioldsphpot=sphpot
                ioldradpot=radpot
                oldradcut=radcut
                radpot=30
                sphpot=302
                radcut=15
            end if
            write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
            write(*,*)
            write(*,*) " Atom index       OW f+          OW f-          OW f0          OW DD"
            call gen1cintgrid(gridatmorg,iradcut)
            sumpos=0
            sumneg=0
            do iatm=1,ncenter
	            gridatm%value=gridatmorg%value !Weight in this grid point
	            gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	            gridatm%y=gridatmorg%y+a(iatm)%y
	            gridatm%z=gridatmorg%z+a(iatm)%z
	            !Calculate orbital weighted f+ and f- at various atomic grids
                OWfposgrid=0
                OWfneggrid=0
	            !$OMP parallel do shared(OWfposgrid,OWfneggrid) private(ipt,imo,wfnval) num_threads(nthreads)
	            do ipt=1,radpot*sphpot
                    call orbderv(1,1,nmo,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,wfnval)
                    do imo=idxLUMO,nmo
                        OWfposgrid(ipt)=OWfposgrid(ipt)+expterm(imo)/denomin_pos*wfnval(imo)**2
                    end do
                    do imo=1,idxHOMO
                        OWfneggrid(ipt)=OWfneggrid(ipt)+expterm(imo)/denomin_neg*wfnval(imo)**2
                    end do
	            end do
	            !$OMP end parallel do
	            !Calculate free atomic density to obtain promolecule density
	            promol=0D0
	            do jatm=1,ncenter
		            !$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
		            do ipt=1,radpot*sphpot
			            tmpdens(ipt)=calcatmdens(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,0)
		            end do
		            !$OMP end parallel do
		            promol=promol+tmpdens
		            if (jatm==iatm) selfdens=tmpdens
	            end do
	            !Now we have needed data in hand, calculate integral
                OW_condfpos=0
                OW_condfneg=0
	            do ipt=1,radpot*sphpot
		            if (promol(ipt)/=0D0) then
                        wei=selfdens(ipt)/promol(ipt)*gridatm(ipt)%value
			            OW_condfpos = OW_condfpos + wei*OWfposgrid(ipt)
			            OW_condfneg = OW_condfneg + wei*OWfneggrid(ipt)
		            end if
	            end do
                OW_condDD=OW_condfpos-OW_condfneg
                OW_condfzero=(OW_condfpos+OW_condfneg)/2
                sumpos=sumpos+OW_condfpos
                sumneg=sumneg+OW_condfneg
                write(*,"(i6,'(',a')',4f15.5)") iatm,a(iatm)%name,OW_condfpos,OW_condfneg,OW_condfzero,OW_condDD
            end do
            write(*,"(/,' Sum of orbital weighted f+',f12.6)") sumpos
            write(*,"(' Sum of orbital weighted f-',f12.6)") sumneg
            if (iautointgrid==1) then
                sphpot=ioldsphpot
                radpot=ioldradpot
                radcut=oldradcut
            end if
        end if
        
    else if (isel==7) then
        if (allocated(cubmat)) deallocate(cubmat)
        aug3Dold=aug3D
        aug3D=3 !Commonly 3 Bohr extension distance is adequate
        call setgrid(1,inouse)
        aug3D=aug3Dold
        allocate(OW_fpos(nx,ny,nz),OW_fneg(nx,ny,nz),cubmat(nx,ny,nz))
        iolduserfunc=iuserfunc
        write(*,*)
        write(*,*) "Calculating grid data of orbital-weighted f+..."
        iuserfunc=95
        call savecubmat(100,0,0)
        OW_fpos=cubmat
        write(*,*)
        write(*,*) "Calculating grid data of orbital-weighted f-..."
        iuserfunc=96
        call savecubmat(100,0,0)
        OW_fneg=cubmat
        iuserfunc=iolduserfunc
        
        sur_value=0.001D0
        do while(.true.)
            write(*,*)
            write(*,"(' -1 Set the value multiplied to all grid data, current:',f12.6)") cubfac
            write(*,*) "0 Return"
            write(*,*) "1 Visualize isosurface of orbital-weighted f+"
            write(*,*) "2 Visualize isosurface of orbital-weighted f-"
            write(*,*) "3 Visualize isosurface of orbital-weighted f0"
            write(*,*) "4 Visualize isosurface of orbital-weighted dual descriptor"
            write(*,*) "5 Export grid data of orbital-weighted f+ as OW_f+.cub"
            write(*,*) "6 Export grid data of orbital-weighted f- as OW_f-.cub"
            write(*,*) "7 Export grid data of orbital-weighted f0 as OW_f0.cub"
            write(*,*) "8 Export grid data of orbital-weighted dual descriptor as OW_DD.cub"
            read(*,*) isel2
            if (isel2==-1) then
                write(*,*) "Input the value, e.g. 0.325"
                read(*,*) cubfac
            else if (isel2==0) then
                deallocate(OW_fpos,OW_fneg)
                exit
            else
                if (isel2==1.or.isel2==5) cubmat=OW_fpos
                if (isel2==2.or.isel2==6) cubmat=OW_fneg
                if (isel2==3.or.isel2==7) cubmat=(OW_fpos+OW_fneg)/2
                if (isel2==4.or.isel2==8) cubmat=OW_fpos-OW_fneg
                cubmat=cubmat*cubfac
                if (isel2<=4) then
		            call drawisosurgui(1)
                else
                    if (isel2==5) open(10,file="OW_f+.cub",status="replace")
                    if (isel2==6) open(10,file="OW_f-.cub",status="replace")
                    if (isel2==7) open(10,file="OW_f0.cub",status="replace")
                    if (isel2==8) open(10,file="OW_DD.cub",status="replace")
                    call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
                    close(10)
                    write(*,*) "Exporting finished!"
                end if
            end if
        end do
        
    else if (isel==8) then !Nucleophilic and electrophilic superdelocalizabilities
		!Also see http://openmopac.net/Manual/super.html. MOPAC can directly calculate by using "super" keyword, its results are given in eV^-1
        call getHOMOidx
        idxLUMO=idxHOMO+1
        alpha_parm=(MOene(idxHOMO)+MOene(idxLUMO))/2
        if (.not.allocated(atmcomp)) allocate(atmcomp(ncenter,nmo))
        call gen_orbatmcomp_space(1,atmcomp,1,nmo,1,1)
        write(*,*)
        write(*,"(' Schuurmann'' alpha parameter:',f16.6,' Hartree')") alpha_parm
        write(*,*) "  All superdelocalizabilities given below are in Hartree^-1"
        write(*,*) "D_N: Nucleophilic superdelocalizability"
        write(*,*) "D_E: Electrophilic superdelocalizability"
        write(*,*) "D_N_0: Nucleophilic superdelocalizability without alpha parameter"
        write(*,*) "D_E_0: Electrophilic superdelocalizability without alpha parameter"
        write(*,*)
        write(*,*) "    Atom           D_N              D_E              D_N_0            D_E_0"
        sumD_N=0
        sumD_E=0
        sumD_N_0=0
        sumD_E_0=0
        do iatm=1,ncenter
			D_E=0
            D_N=0
			D_E_0=0
            D_N_0=0
			do imo=1,idxHOMO
				D_E = D_E + atmcomp(iatm,imo)/(MOene(imo)-alpha_parm)
				D_E_0 = D_E_0 + atmcomp(iatm,imo)/MOene(imo)
            end do
            D_E=D_E*2
            D_E_0=D_E_0*2
			do imo=idxLUMO,nmo
				D_N = D_N + atmcomp(iatm,imo)/(alpha_parm-MOene(imo))
                if (MOene(imo)/=0) then !When using diffuse functions, highest some MOs may be artificially filled and has energy of zero
					D_N_0 = D_N_0 + atmcomp(iatm,imo)/(-MOene(imo))
                end if
            end do
            D_N=D_N*2
            D_N_0=D_N_0*2
            write(*,"(i6,'(',a')',4f17.5)") iatm,a(iatm)%name,D_N,D_E,D_N_0,D_E_0
			sumD_N = sumD_N + D_N
			sumD_E = sumD_E + D_E
			sumD_N_0 = sumD_N_0 + D_N_0
			sumD_E_0 = sumD_E_0 + D_E_0
        end do
        write(*,"(/,' Sum of D_N:  ',f16.5,' /Hartree')") sumD_N
        write(*,"(' Sum of D_E:  ',f16.5,' /Hartree')") sumD_E
        write(*,"(' Sum of D_N_0:',f16.5,' /Hartree')") sumD_N_0
        write(*,"(' Sum of D_E_0:',f16.5,' /Hartree')") sumD_E_0
        
    end if
end do
    
end subroutine