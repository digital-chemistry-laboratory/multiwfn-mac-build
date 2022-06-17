module boysfunc
implicit none
include "boysfunc_data1.h"
contains
    subroutine initboys()
        ! ================================================
include "boysfunc_data2.h"
    end subroutine
    ! Best way to calculate Boys function.
    function boys(x,n)
        real*8 :: boys
        real*8,intent(in) :: x
        integer,intent(in) :: n
        integer :: idx1,idx2,idx3,idx4,idx5
        real*8 :: x1,x2,x3,x4,x5
        real*8 :: y1,y2,y3,y4,y5
        real*8 :: y12,y23,y34,y45
        real*8 :: y123,y234,y345
        real*8 :: y1234,y2345
        ! ================================================
        ! Accurate method.
        if (n > MaxN) then
            boys = accboys(x,n)
        else
            ! Asympotic approximation.
            if (x > AsymX) then
                boys = BoysAsympoticConstants(n+1)*(x**(-n-0.5))
            ! Fast method.
            else
                ! Look up.
                idx1 = x*100+1; x1 = (idx1-1)/100.0; y1 = BoysTable(idx1,n+1)
                idx2 = idx1+1; x2 = x1+dX1; y2 = BoysTable(idx2,n+1)
                idx3 = idx1+2; x3 = x1+dX2; y3 = BoysTable(idx3,n+1)
                idx4 = idx1+3; x4 = x1+dX3; y4 = BoysTable(idx4,n+1)
                idx5 = idx1+4; x5 = x1+dX4; y5 = BoysTable(idx5,n+1)
                ! Neville interpolation.
                y12 = (x-x2)*(y1-y2)/(x1-x2)+y2
                y23 = (x-x3)*(y2-y3)/(x2-x3)+y3
                y34 = (x-x4)*(y3-y4)/(x3-x4)+y4
                y45 = (x-x5)*(y4-y5)/(x4-x5)+y5                
                y123 =(x-x3)*(y12-y23)/(x1-x3)+y23
                y234 =(x-x4)*(y23-y34)/(x2-x4)+y34
                y345 =(x-x5)*(y34-y45)/(x3-x5)+y45
                y1234 = (x-x4)*(y123-y234)/(x1-x4)+y234
                y2345 = (x-x5)*(y234-y345)/(x2-x5)+y345
                boys = (x-x5)*(y1234-y2345)/(x1-x5)+y2345
            endif
        endif
    end function
    ! Calculate Boys function accurately but SLOWLY.
    function accboys(x,n)
        real*8 :: accboys
        real*8,intent(in) :: x
        integer,intent(in) :: n
        integer,parameter :: NIt = 4
        real*8 :: tdts(NIt)
        real*8 :: tresults(NIt)
        integer :: ngs
        integer :: i,n2,j
        real*8 :: dt,t,tresult,tresult12,tresult23,tresult34,tresult123,tresult234
        ! ================================================
        n2 = n*2
        ngs = 25
        ! Naive grid integral.
        do i = 1,NIt
            ngs = ngs*4
            dt = 1./ngs
            if (n2 == 0) then
                tresult = 0.5
            else
                tresult = 0
            endif
            t = dt
            do j = 1,ngs
                tresult = tresult+exp(-x*t*t)*(t**n2)
                t = t+dt
            enddo
            tresult = tresult+exp(-x*t*t)*(t**n2)/2
            tresult = tresult*dt
            tdts(i) = dt
            tresults(i) = tresult
        enddo
        ! Neville extrapolation to zero step.
        tresult12 = (-tdts(2))*(tresults(1)-tresults(2))/(tdts(1)-tdts(2))+tresults(2)
        tresult23 = (-tdts(3))*(tresults(2)-tresults(3))/(tdts(2)-tdts(3))+tresults(3)
        tresult34 = (-tdts(4))*(tresults(3)-tresults(4))/(tdts(3)-tdts(4))+tresults(4)
        tresult123 = (-tdts(3))*(tresult12-tresult23)/(tdts(1)-tdts(3))+tresult23
        tresult234 = (-tdts(4))*(tresult23-tresult34)/(tdts(2)-tdts(4))+tresult34
        accboys = (-tdts(4))*(tresult123-tresult234)/(tdts(1)-tdts(4))+tresult234
    end function
end module
