module mod_2F2
  use util, only: to_lower => struc2lc
  implicit none
  private
  integer, parameter :: GENERAL_EXACT_EVALUATION = 1, &
                        INTERPOLATION_EVALUATION = 2, &
                        SPECIAL_EXACT_EVALUATION = 3
  ! normally, it should be 0 (fractional integrals) or 1 (fractional derivatives)
  integer :: p = 2000
  ! order of fractional derivative
  real(8) :: alpha_
  ! used algorithm
  integer :: algorithm_ = -1
  ! interface to C function; see 2F2.h
  ! note, this routine returns 0 < hyp2F2 <= 1
  ! z should be <= 0
  interface
    real(8) function hyp2F2(a, b, c, d, z) bind(C, name="hyp2F2")
      real(8), intent(in), value :: a, b, c, d, z
    end function hyp2F2
  end interface
  ! internal global arrays for interpolation algorithm
  integer, parameter :: n_2F2 = 101, k_2F2 = 10
  real(8) eprr_2F2(0:n_2F2 - 1)
  real(8) f_2F2(0:n_2F2 - 1, 0:8), Bs_2F2(0:n_2F2 - 1, 0:8)
  real(8) t_2F2(n_2F2 + k_2F2)
  !
  public :: set_alpha_level, GTO_fractional_integral
  public :: p
contains
  subroutine set_alpha_level()
    character(len=80) :: tempstr
    do while (.true.)
      write (*, *) "Type `get` for getting information about current fractional derivative parameters"
      write (*, *) "Type `set_p` for updating fractional derivative parameters with specific integer order of internal derivatives"
      write (*, *) "Type `set_auto` for updating fractional derivative parameters with automatic integer order of internal derivatives"
      read (*, *) tempstr
      call to_lower(tempstr)
      select case(trim(tempstr))
        case("get")
          call get_alpha_level()
          return
        case("set_p")
          write (*, *) "Input fractional order of derivative (alpha); it should be <= 1"
          read (*, *) alpha_
          write (*, *) "Input integer order of internal derivative (p); it can be 0 or 1"
          read (*, *) p
          if (alpha_ > p) then
            write (*, *) "alpha can not be larger than p"
            write (*, *) ""
            cycle
          end if
          if (p /= 0 .and. p /= 1) then
            write (*, *) "p can be only 0 or 1"
            write (*, *) ""
            cycle
          end if
        case("set_auto")
          write (*, *) "Input fractional order of derivative (alpha); it should be <= 1"
          read (*, *) alpha_
          p = ceiling(alpha_)
          if (p < 0) p = 0
          write (*, '(A,I0)') "Automatically chosen integer order of internal derivative (p) is ", p
        case default
          write (*, *) "Try again..."
          cycle
      end select
      exit
    end do
    write (*, *) "Type `exact` for enabling evaluation of general analytic form of fractional derivative"
    write (*, *) "Type `approx` for interpolation of fractional derivative"
    write (*, *) "Type `special` for enabling evaluation of specific analytic form of fractional derivative"
    write (*, *) "For special, only several pairs are available:"
    write (*, *) "    p      alpha"
    write (*, *) "    1       1.0"
    write (*, *) "    0       0.0"
    do while (.true.)
      read (*, *) tempstr
      call to_lower(tempstr)
      select case(trim(tempstr))
        case("exact")
          algorithm_ = GENERAL_EXACT_EVALUATION
        case("approx")
          algorithm_ = INTERPOLATION_EVALUATION
          call prepare_interpolation()
        case("special")
          algorithm_ = SPECIAL_EXACT_EVALUATION
        case default
          write (*, *) "Try again..."
          cycle
      end select
      exit
    end do
    call get_alpha_level()
  end subroutine set_alpha_level
  subroutine get_alpha_level()
    write (*, *) "Fractional derivative parameters:"
    write (*, *) "Alpha: ", alpha_
    write (*, '(A,I0,A,I0)') "Order of derivative is ", p, "; should be ", ceiling(alpha_)
    if (algorithm_ == GENERAL_EXACT_EVALUATION) then
      write (*, *) "Evaluation of fractional derivative using general analytic form"
    else if (algorithm_ == INTERPOLATION_EVALUATION) then
      write (*, *) "Evaluation of fractional derivative using interpolation"
    else if (algorithm_ == SPECIAL_EXACT_EVALUATION) then
      write (*, *) "Evaluation of fractional derivative using specific analytic form"
    else
      write (*, *) "Unknown algorithm"
    end if
    write (*, *)
  end subroutine get_alpha_level
  ! This routine evaluates the following expression:
  ! \frac{1}{\Gamma(p-\alpha)} \int\limits_0^1 (1-t)^{p-\alpha-1} t^n exp(x t^2) dt
  ! using general, approximational or special algorithms
  ! n is integer, principal quantum number
  ! x is real(8) <= 0
  real(8) function GTO_fractional_integral(n, x)
    use bspline_sub_module, only: db1val
    real(8), parameter :: one = 1e0_8, zero = 0e0_8
    integer, intent(in) :: n
    real(8), intent(in) :: x
    real(8) :: top, down
    real(8) :: gtop, gdown
    real(8) :: val2F2
    integer :: iout, inbvx
    if (x > zero) then
      write (*, *) "x is ", x
      error stop "x greater than zero!"
    end if
    val2F2 = -one
    GTO_fractional_integral = 0e0_8
    if (algorithm_ == GENERAL_EXACT_EVALUATION) then
      top = real(n + 1, kind=8)
      down = top + real(p, kind=8) - alpha_
      gtop = gamma(top)
      gdown = gamma(down)
      val2F2 = hyp2F2(top, top + one, down, down + one, x)
      GTO_fractional_integral = gtop / gdown * val2F2
    else if (algorithm_ == INTERPOLATION_EVALUATION) then
      top = real(n + 1, kind=8)
      down = top + real(p, kind=8) - alpha_
      gtop = gamma(top)
      gdown = gamma(down)
      inbvx = 1
      ! normally, x shall be passed, but here some logic is broken since Bspline expects increasing values of x
      call db1val(-x, 0, t_2F2, n_2F2, k_2F2, Bs_2F2(:, n), val2F2, iout, inbvx)
      GTO_fractional_integral = gtop / gdown * val2F2
    else if (algorithm_ == SPECIAL_EXACT_EVALUATION) then
      select case(p)
        case(1)
          if (alpha_ == one) then
            GTO_fractional_integral = exp(x)
          end if
        case(0)
          if (alpha_ == zero) then
            GTO_fractional_integral = exp(x)
          end if
      end select
      if(GTO_fractional_integral /= -one) return
      write (*, *) "alpha: ", alpha_, " p: ", p
      error stop "No specialization for given alpha and p"
    else
      write (*, '(A,I0)') "Algorithm internal value is ", algorithm_
      error stop "Algorithm is not implemented!"
    end if
  end function GTO_fractional_integral
  ! this routine fills arrays for interpolation evaluation of 2F2 function
  subroutine prepare_interpolation()
    use bspline_sub_module, only: db1ink
    real(8), parameter :: one = 1e0_8, zero = 0e0_8
    real(8) :: top, down
    integer :: i, j
    integer :: iout
    do i = 0, n_2F2 - 1
      eprr_2F2(i) = exp(0.1_8 * real(i, 8)) - one
      do j = 0, 8
        if (i > 0 .and. f_2F2(max(i - 1, 0), j) < 1e-40_8) then
          f_2F2(i, j) = zero
        else
          top = real(j + 1, kind=8)
          down = top + real(p, kind=8) - alpha_
          ! we evaluate negative eprr
          f_2F2(i, j) = hyp2F2(top, top + one, down, down + one, -eprr_2F2(i))
        end if
      end do
      write(*,*) i, eprr_2F2(i), f_2F2(i, 0:4)
    end do
    do j = 0, 8
      ! here, eprr_2F2 is from 0 to +inf
      call db1ink(eprr_2F2, n_2F2, f_2F2(:, j), k_2F2, 0, t_2F2, Bs_2F2(:, j), iout)
    end do
  end subroutine prepare_interpolation
end module mod_2F2
