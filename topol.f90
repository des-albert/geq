!-----------------------------------------------------------------------
!     topol.f90
!     Program for identifying region with (g(j,n)-psicon)) > 0
!     connected to the point jaxis,naxis-
!     Outside this region g(j,n) is set equal to 0.0,
!     inside to g(j,n) - psicon.
!     Return code lcod: if > 0, no closed surface existing
!
!     if naxis = 1, top-bottom symmetry is assumed
!
!-----------------------------------------------------------------------
subroutine topol(lcod)
    use solver_mod

    integer(ik) :: lcod, lsym, j, n, jma, l, jmax, jnu, ji, jmin, j1, j2, jmi
    integer(ik), dimension(2) :: nmx
    integer(ik), dimension(2) :: ndi = (/1, -1/)

    
    lsym = 1
    nmx(2) = 0
    if (naxis /= 1) lsym = 2
    lcod = 0

    do j = 1, Mr
        do n = 1, Nz
            g(j,n) = (g(j,n) - psicon)
            if (g(j,n) <= 0.0_rk) g(j,n) = 0.0_rk
        end do
    end do

    n = naxis
    jma = jaxis
    if (g(jma,n) == 0.0_rk) then
        lcod = 5
        print *, 'Plasma disappeared of region, probably squeezed of f at center'
    else
        do l = 1, lsym
            jma = jaxis
            n = naxis
            do j = jma, Mr
                if (g(j,n) == 0.0_rk) goto 10
            end do
            go to 100
10          jmax = j - 1
            jnu = Mr - jmax + 1
            do ji = jnu, Mr
                j = Mr - ji + 1
                if (g(j,n) == 0.0_rk) go to 20
            end do
            go to 110
20          jmin = j + 1
            do while (.true.)
                j1 = jmin - 1
                j2 = jmax + 1
                do j = 1, j1
                    g(j,n) = 0.0_rk
                end do
                do j = j2, Mr
                    g(j,n) = 0.0_rk
                end do
                if (jmax <= jmin) go to 90
                n = n + ndi(l)
                if (.not. ((n < Nz) .and. (n > 1))) go to 140
                if (g(jmax,n) == 0.0_rk) then
                    jnu = Mr - jmax + 1
                    do ji = jnu, Mr
                        j = Mr - ji + 1
                        if (g(j,n) /= 0) go to 30
                    end do
                    go to 70
    30             jmax = j
                else
                    jma = jmax
                    do j = jma, Mr
                        if (g(j,n) == 0.0_rk) go to 40
                    end do
                    go to 120
    40             jmax = j - 1
                end if
                jmi = jmin
                if (g(jmin,n) == 0.0_rk) then
                    do j = jmi, Mr
                        if (g(j,n) /= 0.0_rk) go to 50
                    end do
                    go to 80
    50             jmin = j
                else
                    jnu = Mr - jmin + 1
                    do ji = jnu, Mr
                        j = Mr - ji + 1
                        if (g(j,n) == 0.0_rk) go to 60
                    end do
                    go to 130
    60             jmin = j + 1
                end if
            end do
70          n = n - ndi(l)
            go to 90
80          n = n - 1
90          nmx(l) = n + ndi(l)
        end do

        do n = 1, Nz
            if ((nmx(1) - n)*(nmx(2) - n) >= 0.0_rk) then
                do j = 1, Mr
                    g(j,n) = 0.0_rk
                end do
            end if
        end do

        return
100     print *, 'Plasma runs out of grid on outside at axis'
        lcod = 1
        go to 150
110     print *, 'Plasma runs out of grid on inside at axis'
        lcod = 6
        go to 150

120     print *, 'Plasma runs out of grid on outside for n = ',n
        lcod = 3
        go to 150
130     print *, 'Plasma runs out of grid on inside for n = ',n 
        lcod = 4
        go to 150 
140     print *, 'Plasma runs out on top immersing probably top conductor'
        lcod = 2
              
    endif

150 print *, ' Case abondoned'
    return

end subroutine topol    