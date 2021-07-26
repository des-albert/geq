!-----------------------------------------------------------------------
!
!     plotit.f90
!
!-----------------------------------------------------------------------
subroutine plotit ()

    use solve_mod

    real, dimension(:), allocatable :: rf, rr, rz
    real :: plev

    allocate( rf(MN), rr(Mr), rz(Nz))

! Convert to real for dislin package    

    rf(1:MN) = real(f(1:MN))
    rr(1:Mr) = real(R(1:Mr))
    rz(1:Nz) = real(Z(1:Nz))


!   Flux contours 

    call setpag('ps4l')
    call disini
    call complx
    call pagera
    call pagfll(255)
    call color('blue')
    call axslen(1120,1600)
    call graf(0., 14., 0., 2., -10., 10., -10., 2.)
    call xaxgit
    call color('green')
    call rlrec(1.3135,6.0605,.749,1.979)
    call rlrec(1.3135,4.0325,.749,1.979)
    call rlrec(1.3135,2.0035,.749,1.979)
    call rlrec(1.3135,-0.0245,.749,1.979)
    call rlrec(1.3135,-2.0535,.749,1.979)
    call rlrec(1.3135,-4.0815,.749,1.979)
    call rlrec(3.47,8.045,.968,.976)
    call rlrec(7.9845,6.8275,.649,.595)
    call rlrec(11.581,3.8275,.708,1.125)
    call rlrec(11.5805,-1.6805,.649,1.125)
    call rlrec(7.985,-6.2575,.82,.945)
    call rlrec(3.4705,-7.069,1.633,.976)
    call color('blue');

    do i= 2, 16
        plev = real(psicon + (i-1)*(fabs - psicon)/15.0_rk)
        call contur(rr, Mr, rz, Nz, rf, plev)
    end do
    call color('red')
    call contur(rr, Mr, rz, Nz, rf, real(psicon))
    call endgrf
    call disfin

!   Current and Pressure profile

    call disini
    call complx
    call pagera
    call pagfll(255)
    call color('blue')
    call polcrv('linear')
    call titlin("Current Density & Pressure", 1)
    call graf(3., 10., 3., 1., 0., 2., 0., 1.)
    call title
    call name('Major Radius (m)','X')
    call color('red')
    call curve(rr, real(cjt), Mr)
    call color('green')
    call curve(rr, real(pr), Mr)
    call disfin


end subroutine plotit