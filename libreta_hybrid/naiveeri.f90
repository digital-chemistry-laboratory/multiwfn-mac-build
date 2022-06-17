! J. Chem. Theory Comput. 2018, 14, 572.
! Author: Zhang, Jun

! Calculate (ab|cd) for un-normalized, primitive GTOs using naive Rys algorithm.
subroutine naive_eri( &
    a, RA, lxA, lyA, lzA, &
    b, RB, lxB, lyB, lzB, &
    c, RC, lxC, lyC, lzC, &
    d, RD, lxD, lyD, lzD, &
    xyzG)
use ryspoly
implicit none        
    real*8,intent(in) :: a, RA(3)
    integer,intent(in) :: lxA, lyA, lzA
    real*8,intent(in) :: b, RB(3)
    integer,intent(in) :: lxB, lyB, lzB
    real*8,intent(in) :: c, RC(3)
    integer,intent(in) :: lxC, lyC, lzC
    real*8,intent(in) :: d, RD(3)
    integer,intent(in) :: lxD, lyD, lzD
    real*8,intent(out) :: xyzG
    real*8,parameter :: PI = 3.14159265358979323846264338327950
    real*8 :: p, mu, RAB(3), RP(3), RPA(3), RPB(3), KAB
    real*8 :: q, nu, RCD(3), RQ(3), RQC(3), RQD(3), KCD
    real*8 :: alpha, RPQ(3), norm, normG
    integer :: gamma, i
    real*8 :: beta, roots(25), weights(25)
    real*8 :: xG, yG, zG
    ! ================================================
    ! GTOs and their overlaps.
    p = a+b
    mu = a*b/p
    RAB = RA-RB
    RP = (RA*a+RB*b)/p
    RPA = RP-RA
    RPB = RP-RB
    KAB = exp(-mu*dot_product(RAB, RAB))
    q = c+d
    nu = c*d/q
    RCD = RC-RD
    RQ = (RC*c+RD*d)/q
    RQC = RQ-RC
    RQD = RQ-RD
    KCD = exp(-nu*dot_product(RCD, RCD))
    RPQ = RP-RQ
    alpha = p*q/(p+q)
    norm = 1.
    normG = 2*(PI**2.5)/p/q/sqrt(p+q)
    ! Calculate Rys roots and weights.
    gamma = (lxA+lxB+lxC+lxD+lyA+lyB+lyC+lyD+lzA+lzB+lzC+lzD)/2+1
    beta = alpha*dot_product(RPQ, RPQ)
    if (gamma <= 3) call calcryspoly_123(gamma, beta, roots, weights) 
    if (gamma == 4) call calcryspoly_4(beta, roots, weights) 
    if (gamma == 5) call calcryspoly_5(beta, roots, weights) 
    if (gamma >= 6) call calcryspoly_ge6(gamma, beta, roots, weights) 
    ! Loop over roots and weights.
    xyzG = 0
    do i = 1, gamma, 1
        call naive_eri_core(lxA, lxB, lxC, lxD, RPA(1), RPB(1), RQC(1), RQD(1), p, q, RPQ(1), alpha, roots(i), xG)
        call naive_eri_core(lyA, lyB, lyC, lyD, RPA(2), RPB(2), RQC(2), RQD(2), p, q, RPQ(2), alpha, roots(i), yG)
        call naive_eri_core(lzA, lzB, lzC, lzD, RPA(3), RPB(3), RQC(3), RQD(3), p, q, RPQ(3), alpha, roots(i), zG)
        xyzG = xyzG+xG*yG*zG*weights(i)
    enddo
    xyzG = xyzG*normG*KAB*KCD
end subroutine

subroutine naive_eri_core(lA, lB, lC, lD, RPA, RPB, RQC, RQD, p, q, RPQ, alpha, root, result)
implicit none
    integer, intent(in) :: lA, lB, lC, lD
    real*8, intent(in) :: RPA, RPB, RQC, RQD, p, q, RPQ, alpha, root
    real*8, intent(out) :: result
    integer lABCD, iABCD, iA, iB, iC, iD
    real*8 :: GABCD(25, 25, 25, 25), CA, CB, CC, CD, CAB, CCD, C0
    ! ================================================
    lABCD = lA+lB+lC+lD    
    CA = RPA-root*RPQ*alpha/p
    CB = RPB-root*RPQ*alpha/p
    CC = RQC+root*RPQ*alpha/q
    CD = RQD+root*RPQ*alpha/q
    CAB = 0.5/p*(1.-alpha/p*root)
    CCD = 0.5/q*(1.-alpha/q*root)
    C0 = 0.5/(p+q)*root
    GABCD(1, 1, 1, 1) = 1
    do iABCD = 1, lABCD, 1
        do iA = 0, lA, 1
            do iB = 0, lB, 1
                do iC = 0, lC, 1
                    iD = iABCD-iA-iB-iC
                    if (iD < 0) cycle
                    if (iA >= 1) then
                        GABCD(iA+1, iB+1, iC+1, iD+1) = CA*GABCD(iA, iB+1, iC+1, iD+1)
                        if (iA >= 2) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                            GABCD(iA+1, iB+1, iC+1, iD+1)+CAB*(iA-1)*GABCD(iA-1, iB+1, iC+1, iD+1)
                        if (iB >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                            GABCD(iA+1, iB+1, iC+1, iD+1)+CAB*iB*GABCD(iA, iB, iC+1, iD+1)
                        if (iC >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                            GABCD(iA+1, iB+1, iC+1, iD+1)+C0*iC*GABCD(iA, iB+1, iC, iD+1)
                        if (iD >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                            GABCD(iA+1, iB+1, iC+1, iD+1)+C0*iD*GABCD(iA, iB+1, iC+1, iD)
                    else
                        if (iB >= 1) then
                            GABCD(iA+1, iB+1, iC+1, iD+1) = CB*GABCD(iA+1, iB, iC+1, iD+1)
                            if (iB >= 2) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                GABCD(iA+1, iB+1, iC+1, iD+1)+CAB*(iB-1)*GABCD(iA+1, iB-1, iC+1, iD+1)
                            if (iA >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                GABCD(iA+1, iB+1, iC+1, iD+1)+CAB*iA*GABCD(iA, iB, iC+1, iD+1)
                            if (iC >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                GABCD(iA+1, iB+1, iC+1, iD+1)+C0*iC*GABCD(iA+1, iB, iC, iD+1)
                            if (iD >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                GABCD(iA+1, iB+1, iC+1, iD+1)+C0*iD*GABCD(iA+1, iB, iC+1, iD)
                        else
                            if (iC >= 1) then
                                GABCD(iA+1, iB+1, iC+1, iD+1) = CC*GABCD(iA+1, iB+1, iC, iD+1)
                                if (iC >= 2) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                    GABCD(iA+1, iB+1, iC+1, iD+1)+CCD*(iC-1)*GABCD(iA+1, iB+1, iC-1, iD+1)
                                if (iA >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                    GABCD(iA+1, iB+1, iC+1, iD+1)+C0*iA*GABCD(iA, iB+1, iC, iD+1)
                                if (iB >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                    GABCD(iA+1, iB+1, iC+1, iD+1)+C0*iB*GABCD(iA+1, iB, iC, iD+1)
                                if (iD >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                    GABCD(iA+1, iB+1, iC+1, iD+1)+CCD*iD*GABCD(iA+1, iB+1, iC, iD)
                            else
                                GABCD(iA+1, iB+1, iC+1, iD+1) = CD*GABCD(iA+1, iB+1, iC+1, iD)
                                if (iD >= 2) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                    GABCD(iA+1, iB+1, iC+1, iD+1)+CCD*(iD-1)*GABCD(iA+1, iB+1, iC+1, iD-1)
                                if (iA >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                    GABCD(iA+1, iB+1, iC+1, iD+1)+C0*iA*GABCD(iA, iB+1, iC+1, iD)
                                if (iB >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                    GABCD(iA+1, iB+1, iC+1, iD+1)+C0*iB*GABCD(iA+1, iB, iC+1, iD)
                                if (iC >= 1) GABCD(iA+1, iB+1, iC+1, iD+1) = &
                                    GABCD(iA+1, iB+1, iC+1, iD+1)+CCD*iC*GABCD(iA+1, iB+1, iC, iD)
                            endif
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo
    result = GABCD(lA+1, lB+1, lC+1, lD+1)
end subroutine

