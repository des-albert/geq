!----------------------------------------------------------------------
!
!  Initial guess at plasma current
! 
!----------------------------------------------------------------------
subroutine startt()
    use solve_mod
    real(rk) :: cp
    integer :: i, j, nmin, njmax

    cp = 3.0_rk*(totcurr)/(20.0_rk*dr*dz)/abs(real(ndes - naxis,rk))  

    write(*,'(a,f10.4)') 'Start current : ', totcurr

    do i = 1, Mr
        do j = 1, Nz
        g(i,j) = 0.0_rk
        end do
    end do

    nmin = 2*naxis - ndes
    njmax = ndes
    if (njmax .le. nmin) then
        njmax = nmin
        nmin = ndes
    end if
    do j = nmin, njmax
        do i = 1, 5
            g(jaxis-3 + i,j) = (cp*(1.0_rk-(real(j-naxis,rk)/real(ndes - naxis,rk))**2))*dz**2
        end do
    end do
    
    return
end subroutine startt