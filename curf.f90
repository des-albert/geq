!-----------------------------------------------------------------------
!     curf.f90
!
!
!
!  Subroutine to find plasma current distribution
!  using psi distribution stored in functions pp and ff.
!

subroutine curf()

    use solve_mod

    integer(ik) :: j, k, n, mcont, ipoi, lcod
    real(rk), dimension(10) :: s, beh, dint
    real(rk) :: rax, zax, xalp, xalp1, xalp2, xalpr2, xalpz2, rzy, gn, cjm
    real(rk) :: pp, ff, ppin, ffin

    if (mprfg.ne.0) write(*,'(a,f12.5)') 'psi saddle point                =  ',psicon

    call topol(lcod)
!
!
    do k = 1, 10
        dint(k) = 0.0_rk
    enddo

    fmaxa = 0.0_rk
    mcont = 0
    do j = 1, Mr
        do n = 1, Nz
            if (g(j,n) .gt. fmaxa) then
                fmaxa = g(j,n)
                jaxis = j
                naxis = n
                rax = R(jaxis)
                zax = Z(naxis)
            end if
        end do
    end do

    fabs = fmaxa + psicon
    if (mprfg.ne.0) write (*,'(a,f12.5,a,f12.5,a,f12.5)') 'Magnetic Axis radius =', rax, ' height = ', zax, ' psi = ', fabs
    ipoi = 0
    idol = 0
    xalp = 0.0_rk

    do j = 1, Mr
        rzy = R(j)
        xalp1 = 0.0_rk
        do k = 1, 10
            beh(k) = 0.0_rk
        end do
        do n = 1, Nz
            xalp2 = 0.0_rk
            xalpr2 = 0.0_rk
            xalpz2 = 0.0_rk
            if (g(j,n) .gt. 0.0_rk) then
                if (((j .gt. 1) .and. (j .lt. Mr)) .and. ((n .gt. 1) .and. (n .lt. Nz))) then
                    xalpr2 = ((g(j+1, n) - g(j-1, n))/(2.0_rk*dr))**2
                    xalpz2 = ((g(j, n+1) - g(j, n-1))/(2.0_rk*dz))**2

                    if (g(j-1, n) .le. 0.0_rk)then
                        xalpr2 = ((g(j+1, n) - g(j,n))/dr)**2
                    end if
                    if (g(j+1, n) .le. 0.0_rk) then
                        xalpr2 = ((g(j,n) - g(j-1,n))/dr)**2
                    end if
                    if (g(j, n-1) .le. 0.0_rk) then
                        xalpz2 = ((g(j,n+1) - g(j,n))/dz)**2
                    end if
                    if (g(j, n+1) .le. 0.0_rk) then
                        xalpz2 = ((g(j,n) - g(j,n-1))/dz)**2
                    end if

                    xalp2 = xalp2 + (xalpr2 + xalpz2)/rzy
                    ipoi = ipoi + 1
                else
                    idol = idol + 1
                end if

                gn = g(j,n)/fmaxa
                s(1) = pp(gn, alfac)*rzy
                s(2) = ff(gn, alfac)/rzy
                s(3) = 1.0_rk/rzy**2
                s(4) = ppin(gn, alfac)*fmaxa
                s(5) = ffin(gn, alfac)/rzy**2*fmaxa
                s(6) = 1.0_rk
                s(7) = 2.0_rk*pi*rzy*ppin(gn, alfac)*fmaxa
                s(8) = 2.0_rk*pi*rzy
                s(9) = 0.0_rk
                s(10) = 0.0_rk
                do k = 1, 10
                    beh(k) = beh(k) + s(k) 
                end do
                
            end if
        end do
    
        xalp = xalp + xalp1
        do k = 1, 10
            dint(k) = dint(k) + beh(k)
        end do                 
    end do
    do k = 1, 10
        dint(k) = dint(k)*dz*dr
    enddo    
    xalp = 2.0_rk*pi*xalp*dr*dz
    xind = 2.0_rk*xalp/(totcurr*totcurr*Rmpl)
    cjm = 0.0_rk
    c0pr = betapol*totcurr**2/(8.0_rk*pi*dint(4))
    c0btr = (totcurr - c0pr*dint(1))/dint(2)

    if(mprfg.ne.0) then
        write(*,'(a,e12.5,a,e12.5)') 'c0pr = ', c0pr, ' c0btr = ', c0btr 
        write(*,'(a,e12.5,a,e12.5)')' Plasma Area = ',dint(6), ' Plasma Volume = ', dint(8)
    end if

    do j = 1, Mr
        gn = g(j,naxis)/fmaxa
        pr(j) = c0pr*ppin(gn, alfac)*fmaxa
        bt2(j) = c0btr*ffin(gn, alfac)/R(j)**2*fmaxa
    end do

    do j = 1, Mr
        do n = 1, Nz
            gn = g(j,n)/fmaxa
            g(j,n) = (c0pr*pp(gn, alfac)*R(j)**2 + c0btr*ff(gn, alfac))*dz**2
        end do
    end do

    pint = c0pr*dint(4)
    pintvo = c0pr*dint(7)
    bt2int = c0btr*dint(5)
    betap = 2.0_rk*pint/(totcurr**2*dint(6))
    betapvo = 2.0_rk*pintvo/(totcurr**2*dint(8))
    betat = dint(3)

    do j = 1, Mr
        cjt(j) = g(j,naxis)/r(j)
        if (abs(cjt(j)) .gt. cjm) cjm = abs(cjt(j))
    end do
    if (cjm .lt. dz**2) cjm = 1.0_rk
    do j = 1, Mr
        cjt(j) = cjt(j)/cjm
    end do

    return

end subroutine curf

function pp (x, alfa)

    use solve_mod

    real(rk) :: pp, ff, ppin, ffin, x, alfa

    pp = x**alfa*(2.0_rk - x**alfa)
    return

    entry  ff(x,alfa)
    ff = x**alfa*(2.0_rk - x**alfa)
    return

    entry  ppin(x,alfa)
    ppin = (2.0_rk/(alfa + 1.0_rk))*x**(alfa+1.0_rk)-(1.0_rk/(2.0_rk*alfa+1.0_rk))*x**(2.0_rk*alfa + 1.0_rk)
    return

    entry  ffin(x,alfa)
    ffin = (2.0_rk/(alfa + 1.0_rk))*x**(alfa+1.0_rk)-(1.0_rk/(2.0_rk*alfa+1.0_rk))*x**(2.0_rk*alfa + 1.0_rk)
    return
    
end function pp