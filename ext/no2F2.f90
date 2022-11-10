module no2F2_WIN
contains
  real(8) function hyp2F2(a, b, c, d, z) bind(C, name="hyp2F2")
    real(8), intent(in), value :: a, b, c, d, z
    write(7,*) "Fractional derivatives are not supported!"
    write(7,*) "Please, recompile Multiwfn with its support!"
    error stop
  end function hyp2F2
end module no2F2_WIN