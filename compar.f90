!
!---------------------------------------------------------------------
!
!
!     compar.f90
!
!---------------------------------------------------------------------
!
!
subroutine compar ()

    use solve_mod

    integer(ik) :: j
    real(rk) :: rel, tot, dev, ren

    if (icycle .ge. 1) then
        rel = 0._rk
        do j = 1, Mr
            tot = 0.5_rk*abs(g(j,naxis)) + 0.5_rk*abs(com(j))
            dev = 0.5_rk*abs(g(j,naxis) - com(j))
            ren = (dev/tot)
            if (ren .gt. rel) rel = ren
        end do
        write (*,'(a,e12.5)') ' Relative Error = ', rel
        if (rel .le. error) then
            idecis = 1
            return
        end if
    end if

    do j = 1, Mr
        com(j) = g(j,naxis)
    end do
    idecis = 0
    return

end subroutine compar