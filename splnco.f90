!-----------------------------------------------------------------------
!
!     splnco.f90
!
!     Author : Des Albert
!
!   Purpose:
!               Quadratic spline coefficients for array psi.
!               Algorithm: O. Buneman, J. Comp. Phys. Vol 11, P 250, 1973
!               this routine initializes psi for later calls of intpol.
!
!   Input:     Values to be interpolated in array psi(m1,n1).
!              m1-1 and n1-1 must both be powers of 2 and >=4
!              or divisible by 32.
!
!-----------------------------------------------------------------------

subroutine splnco (psi)
    use solver_mod

    real(rk),dimension(8) ::  a
    real(rk),dimension(*) :: psi
    real(rk) :: a1

    integer(ik) :: i, j, k, l, jh, il, jh1, jh2, lend, ls, mode, init, iend 

    a1 =  - 4.0_rk
    do i = 1, 4
        a(i) = 1.0_rk/a1
        a1 = a1 - 2.0_rk
        a(9-i) = 1.0_rk/a1
        a1 = a1*a1
    end do

!
!   Control for inner loop on columns
!
    jh1 = 1
    lend = Mr*Nm1 + 1
    ls = Mr
    jh2 = min0(16, Mm1/2)
    il = Mm1 - 1


!
!   Loop through columns and rows
! 
20  k = 1
    jh = jh1
    mode = 2
30  j = 2*jh    

    do l = 1, lend, ls
        init = l + jh*mode
        iend = l + il
        do i = init, iend, j
            psi(i) = psi(i) + (psi(i+jh)+psi(i-jh) - 2.0_rk*psi(i))*a(k)
        end do
    end do
    k = k + 1

    select case (mode)
        case(2)
            jh = 2*jh
            if (jh < jh2) go to 30
            mode = 1
            if (k == 5) then
                jh = jh/2
                if (jh >= jh1) go to 30
                if (jh1 == Mr) go to 50
                jh1 = Mr
                lend = Mr
                ls = 1
                jh2 = min0(Mr*16, Mr*(Nm1/2))
                il = (Nm1-1)*Mr
                go to 20
            end if
            k = 9 - k
            go to 30
        case(1) 
            jh = jh/2
            if (jh >= jh1) go to 30
            if (jh1 == Mr) go to 50
!
!        Control for inner loop on rows
!
            jh1 = Mr
            lend = Mr
            ls = 1
            jh2 = min0(Mr*16, Mr*(Nm1/2))
            il = (Nm1 - 1)*Mr
            go to 20
    end select 
50  continue
    return

end subroutine splnco