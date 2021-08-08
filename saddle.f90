!-----------------------------------------------------------------------!
!     saddle.f90
!
!       Purpose:
!               Finds the saddle-point with largest psi value within a rectangle
!   Input:
!       array   q(m,n): values of psi
!       m:      number of points in r direction >=3
!       n:      number of points in z direction >=3
!   Output:
!       irsp:   >=2:    (irsp,izsp,psicon ) r,z index and psi value of saddle-point
!       irsp:   = 0     no saddle-point found
!       irsp:   = -1    m or n < 3
!-------------------------------------------------------------------------------
subroutine saddle()

    use solver_mod

    logical :: sg, gg
    integer(ik) :: i, k, ko, io, l
    real(rk) :: hk, zfrr, zfzz, zfrz

!
!   Seach internal points
!
    irsp = 0
    do i = 2, Nm1
        do k = 2, Mm1
            hk = g(k,i)
            if (.not. (irsp /= 0 .and. psicon >= hk)) then
                zfrr = g(k+1,i) - 2.0_rk*hk + g(k-1,i)
                zfzz = g(k,i+1) - 2.0_rk*hk + g(k,i-1)
                zfrz = g(k+1,i+1) + g(k-1,i-1) - g(k+1,i-1) - g(k-1,i+1)
                if (16._rk*zfrr*zfzz - zfrz**2 < 0.0_rk) then
                    sg = .false.
                    gg = .false.
                    ko = 1
                    io =  - 2
                    do l = 1, 4
                        if (l < 4) io = io + 1
                        if (l == 4) ko = 0
                        if (f1(g(k+ko, i+io), g(k-ko,i-io)) > 0.0_rk) then
                            if (g(k+ko, i+io) /= hk) then
                                if (g(k+ko, i+io) > hk) then
                                    sg = .true.
                                else
                                    gg = .true.
                                end if
                                if (sg .and. gg) go to 10
                            endif
                        end if
                    end do
                    go to 20
10                  irsp = k
                    izsp = i
                    psicon = hk
                end if
            end if
20      continue
        end do
    end do

    return
      
    contains
    
function f1(a,b)
    real(rk) :: f1, a, b
    f1 = (a - hk)*(b - hk)
end function f1

end subroutine saddle