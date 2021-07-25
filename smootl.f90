!-----------------------------------------------------------------------
!
subroutine smootl (y, imax, b, c)

!---------------------------------------------------------------------
    use solve_mod

    real(rk):: b, c
    real(rk), dimension(5) :: sa, ssa 

    integer(ik) :: i, j, imax 

    do i = 1,3
        sa(i) = 0.0_rk
        ssa(i) = 0.0_rk
    end do

    sa(1) = 1.0_rk
    ssa(1) = y(1)
    do  j = 2, imax
        do  i = 1,3
            sa(i) = sa(i) + (j-1)**(i-1)
            ssa(i) = ssa(i) + y(j)*(j-1)**(i-1)
        end do
    end do
    c = (sa(3)*ssa(1) - sa(2)*ssa(2))/(sa(3)*sa(1) - sa(2)**2)
    b = (ssa(1) - c*sa(1))/sa(2)
    
    return
    
    end subroutine smootl