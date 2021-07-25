!-------------------------------------------------------------------------------
!
!     gelg.f90
!
!   Purpose:
!
!       Solves a general system of simultaneous linear equations
!
!   Parameters:
!
!   Input:  
!       r(m,n), m x n matrix of right hand sides ( destroyed )
!               on return, r contains the solution of the equations.
!       a:      the m x m coefficient matrix ( destroyed )
!       m:      the number of equations in the system
!       n:      the number of right hand side vectors   
!       eps:    an input constant which is used as relative tolerance for
!               test on loss of significance
!       ier:    resulting error parameter coded as follows
!               ier  =  0  - no error,
!               ier  =  -1 - no result because of m less than 1 or
!                       pivot element at any elimination step equal to 0,
!               ier  =  k  - warning due to possible loss of significance
!                   indicated at elimination step k+1 where pivot
!                   element was less than or equal to the internal
!                   tolerance eps times absolutely greatest     element
!                   of matrix a.
!   
!       Remarks:
!           Input matrices rr and a are assumed to be stored columnwise
!           in m*n respectively m*m successive storage locations.
!           On return solution matrix r is stored columnwise too.
!           The procedure gives results if the number of equations m is
!           greater than 0 and pivot elements at all elimination steps
!           are different from 0. However warning ier = k - if given
!           indicates possible loss of significance. In case of a well
!           scaled matrix a and appropriate tolerance eps, 
!           ier = k may be interpreted that matrix a has the rank k. 
!           No warning is given in case m = 1.
! 
! 
!  Method:
!           Solution is done by means of gauss-elimination with complete pivoting.
!------------------------------------------------------------------------------- 
subroutine gelg (rr, a, m, n, eps, ier) 
    
    use solve_mod

    integer(ik) m, n, ier, mm, nm, lst, i, j, k, ll, lend, ii, ist
    real(rk), dimension(*) :: a, rr
    real(rk) :: eps, piv, pivi, tb, tol

    if (m .gt. 0) then
!
!     Search for greatest element in matrix a
!
        ier = 0
        piv = 0._rk
        mm = m*m
        nm = n*m
        do l = 1, mm
            tb = abs(a(l))
            if (tb .gt. piv) then
                piv = tb
                i = l
            endif
        enddo
        tol = eps*piv
!
!     a(i) is pivot element. piv contains the absolute value of a(i).
!
!
!     Start elimination loop

        lst = 1
        do k = 1, m
!
!     Test on singularity
!
            if (piv .le. 0) go to 230
            if (ier .eq. 0) then
                if (piv .le. tol) ier = k - 1
            end if
            pivi = 1.0_rk/a(i)
            j = (i - 1)/m
            i = i - j*m - k
            j = j + 1 - k
!
!     i+k is row-index, j+k column-index of pivot element
!     pivot row reduction and row interchange in right hand side r
!
            do l = k, nm, m
                ll = l + i
                tb = pivi*rr(ll)
                rr(ll) = rr(l)
                rr(l) = tb
            enddo           
!
!     Is elimination terminated
!
            if (k .ge. m) go to 180
!
!     Column interchange in matrix a
!
            lend = lst + m - k
            if (j .gt. 0) then
                ii = j*m
                do l = lst, lend
                    tb = a(l)
                    ll = l + ii
                    a(l) = a(ll)
                    a(ll) = tb
                end do
            end if  
!
!     Row interchange and pivot row reduction in matrix a
!
            do l = lst, mm, m
                ll = l + i
                tb = pivi*a(ll)
                a(ll) = a(l)
                a(l) = tb
            end do
!
!     Save column interchange information
!
            a(lst) = j
!
!     Element reduction and next pivot search
! 
            piv = 0.0_rk
            lst = lst + 1
            j = 0
            do ii = lst, lend
                pivi =  - a(ii)
                ist = ii + m
                j = j + 1
                do l = ist, mm, m
                    ll = l - j
                    a(l) = a(l) + pivi*a(ll)
                    tb = abs(a(l))
                    if (tb .gt. piv) then
                        piv = tb
                        i = l
                    end if
                end do
                do l = k, nm, m
                    ll = l + j
                    rr(ll) = rr(ll) + pivi*rr(l)
                end do
            end do
            lst = lst + m
        end do
!
!     End of elimination loop
!
!     Back substitution and back interchange
!
180     if (m .ge. 1) then
            if (m .ne. 1) then
                ist = mm + m
                lst = m + 1
                do i = 2, m
                    ii = lst - i
                    ist = ist - lst
                    l = ist - m
                    l = int(a(l) + 0.5_rk)
                    do j = ii, nm, m
                        tb = rr(j)
                        ll = j
                        do k = ist, mm, m
                            ll = ll + 1
                            tb = tb - a(k)*rr(ll)
                        end do
                        k = j + l
                        rr(j) = rr(k)
                        rr(k) = tb
                    end do
                end do
            end if
        
            return

        end if
    end if

230 ier =  - 1
    return

        
end subroutine gelg
