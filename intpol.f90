!-----------------------------------------------------------------------
!
subroutine intpol (v, x, y, fi, fx, fy, ier)
    use solver_mod

    real(rk), dimension(*):: v
    real(rk),dimension (3) :: af,adf
    real(rk) :: x, y, fi, fx, fy, dx, dy, a1, a2
    integer(ik) :: ier, m, n, mni, ki

    dx = x - int(x + 0.5_rk)
    dy = y - int(y + 0.5_rk) 
    m = int(x - 0.5_rk) 
    n = int(y + 0.5_rk) 
    mni = m + n*Mr
    ier = 1

    if (.not. (m .lt. 0 .or. m .ge. Mr-2 .or. n .lt. 0 .or. n .ge. Nz -1)) then
        ier = 0
        if (n .eq. 0) then
            do i = 1, 3
                ki = mni + 1
                a1 = 2.0_rk*(v(k + Mr) - v(ki))
                a2 = 0.0_rk
                af(i) = p(dy) + v(ki)
                adf(i) = dp(dy)
            end do
        else
            do i = 1, 3
                ki = mni + i
                a1 = v(ki + Mr) + v(ki - Mr) - 2.0_rk*v(ki)
                a2 = v(ki + Mr) - v(ki - Mr)
                af(i) = p(dy) + v(ki)
                adf(i) = dp(dy)
            end do
        end if
        a1 = af(1) + af(3) - 2.0_rk*af(2)
        a2 = af(3) - af(1)
        fi = p(dx) + af(2)
        fx = dp(dx)
        a1 = adf(1) + adf(3) - 2.0_rk*adf(2)
        a2 = adf(3) - adf(1)
        fy = p(dx) + adf(2)
    end if  

    return

    contains

    function p(x)
        real(rk) x,p
        p = 0.5_rk*(a1*x + a2)*x + 0.125_rk*a1
        end function

        function dp(x)
        real(rk) x,dp
        dp = a1*x + 0.5_rk*a2
    end function

end subroutine intpol


subroutine condit (v, rx, zx, itp, b)
    
    use solver_mod

    real(rk) :: rx, zx, x, y, b
    real(rk),dimension(*):: v
    real(rk), dimension(3) :: var
    integer(ik) :: itp, ier

    x = (rx - Rmin)/dr
    y = (zx - Zmin)/dz

    call intpol(v, x, y, var(1), var(3), var(2), ier)
    b = var(itp)

    return

end subroutine condit