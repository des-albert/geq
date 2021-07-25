!----------------------------------------------------------------------
!     flux.f90
!
!   Purpose:
!
!       Solves the Axisymmetric Flux Equation in rectangular Torus
!       cross-section with coordinates r,z:
!       
!       delstar(psi) = -f(r,z)
!
!       with arbitrary mesh ratio dz/dr and prescribed boundary values 
!       for psi,  without symmetry with respect to z=0.
!       the operator delstar is defined as
!
!               delstar(psi)=r*d/dr(1/r*dpsi/dr)+d/dz(dpsi/dz)
!      
!       and is approximated by direct differencing.
!
!       psi:    must contain f(r,z)*dz**2 in the interior and specified
!               boundary values on the sides, first index counting
!               r-direction.
!
!
!
!   Method:  This is an adaption of Buneman's Poisson Solver,S
!      suipr report no. 294, 1969, Stanford University
!----------------------------------------------------------------------

subroutine flux(psi)
    use solve_mod

    real(rk), dimension(:), allocatable :: tc, a
    real(rk), dimension(*) :: psi
    integer(ik) :: nn, itc, lo, l, iu, j1, ju, is, mode, iphase, jd, jh, jt, ji
    integer(ik) :: j, jdm, jhm, jtm, jiu
    real(rk) :: d,b

    allocate(tc(llp), a(llp))

    nn = Nm1
    iu = Mm1 - 1
    itc = 2*Mm1
    tc(itc + nn/2) = 0.0_rk
    lo = nn/2
10  l = lo/2
    tc(itc + l) = sqrt(2.0_rk + tc(itc+lo))
    lo = l
20  tc(itc + nn - l) = -tc(itc+l)
    l = l + 2*lo

    select case ((2*l/nn)*(2*lo-3))

    case (1:)
        go to 10
    case (0) 
        tc(itc + l) = (tc(itc + l + lo) + tc(itc + l - lo))/tc(itc + lo)
        go to 20
    case (:-1)
        do i = 1, iu
            d = alpha/(float(Mm1) + alpha*float(2*i - Mm1))
            tc(i) = ss/(1.0_rk - d)
            tc(i + Mm1) = ss/(1.0_rk + d)
        end do
    end select
    
    lo = nn/2
    j1 = 1 + Mr*(Nm1/nn)
    ju = Nm1*Mr

    do i = j1, ju, Mr
        psi(i + 1) = psi(i + 1) + tc(1)*psi(i)
        psi(iu + i) = psi(iu + i) + tc(iu + Mm1)*psi(i + Mm1)
    end do

    a(Mm1) = 0.0_rk
    a(Mm1 + Mm1) = 0.0_rk
    mode = 2
    is = -1
80  li = 2*lo

    iphase = 2*mode - li/nn
    jd = Mr*nn/li
    jh = jd/2
    jt = jd + jh
    ji = 2*jd
    jo = jd*mode*((1-is)/2) + 1
    do j = jo, ju, ji
        j1 = j + 1
        jdm = jd*is
        jhm = jh*is
        jtm = jt*is
        jiu = j + iu

        select case (iphase)
            case (1)
                do i = j1, jiu
                    a(i-j) = psi(i) + psi(i + jd) + psi(i + jdm)
                    psi(i) = 0.0_rk
                end do
            case (2)
                do i = j1, jiu
                    a(i-j) = psi(i) + psi(i + jd) + psi(i + jdm)
                    psi(i) = 0.5_rk*(psi(i) - psi(i + jh) - psi(i + jhm))
                end do
            case (3)
                do i = j1, jiu
                    a(i-j) = 2.0_rk*psi(i)
                    psi(i) = psi(i + jd) + psi(i + jdm)
                end do
            case (4)
                do i = j1, jiu
                    d = psi(i) - psi(i + jt) - psi(i + jtm)
                    psi(i) = psi(i) - psi(i + jh) - psi(i + jhm) + psi(i + jd) + psi(i + jdm)
                    a(i-j) = d + psi(i)
            end do
        end select

        do l = lo, nn, li
            d = 2.0_rk - tc(itc + l)
            do i = 1, iu
                k = Mm1 - i
                b = 1.0_rk/(d + tc(k) + tc(k + Mm1)*(1.0_rk - a(k + Mr)))
                a(k + Mm1) = tc(k)*b
                a(k) = (a(k) + tc(k + Mm1)*a(k + 1))*b
            end do
            do i = 2, iu
            a(i) = a(i + Mm1)*a(i - 1) + a(i)
            end do
        end do
        do i = j1, jiu
            psi(i) = psi(i) + a(i - j)
        end do
        is = -1
    end do

    select case (iphase)
        case (1)
            return
        case (2)
            lo = 2*lo
            if (lo.lt.nn) go to 80
        case (3:4)
            lo = lo/2
            if (lo.eq.1) mode = 1
        go to 80
    end select

    return
	  
end subroutine flux