!-----------------------------------------------------------------------!
!     xcur.f90
!
!     Stored in cl(i,k) the mutual inductances,
!     in cl(mmax+1,i) the algebraic sum of the winding number of the i-th conductor package
!
!
!-----------------------------------------------------------------------
subroutine xcur(a)

    use solver_mod

    integer :: mmaxp1, mmaxp2, i, j, ll, k, ier, ki
    real(rk), dimension(mpnmax,1) :: a
    real(rk) :: energy

    real(rk), dimension(3) :: bv = (/ -1._rk, 0.0_rk, 0._rk /)

    mmaxp1 = Mmax + 1
    mmaxp2 = Mmax + 2

    do i = 1, Mmax
        do j = i, Mmax
            a(i,j) = cl(i,j)
        end do
        a(i,Mmax+1) = 0.0_rk
        do j = 1, Nmax
            a(i,mmaxp1+j) = bb(j,i)
        end do
        if (icops >= 2) then
            if (icops > 2) then
                a(i,mpnmax) = 0.0_rk
            else
                a(i,mpnmax) = cl(mmaxp1,i)
            end if
        end if
    end do

    a(mmaxp1,mmaxp1) = 0._rk
    if (Nmax > 0) then
        do j = 1, Nmax
            a(mmaxp1, mmaxp1 + j) = bv(ityp(j))
        end do
    end if
    if (icops >= 2) then
        if (icops > 2) then
            a(mmaxp1,mpnmax) = 1.0_rk
        else
            a(mmaxp1,mpnmax) = 0.0_rk
        end if
    end if
    do i = mmaxp2, mpnmax
        do j = i, mpnmax
            a(i,j) = 0.0_rk
        end do
    end do
    do i = 1, mmaxp1
        fk(i) = 0.0_rk
    end do

    if (llmax > 0) then
        do ll = 1, llmax
            do i = 1, Mmax
                do j = i, Mmax
                    a(i,j) = a(i,j) + 2.0_rk*alph*eb(ll,i)*eb(ll,j)
                end do
                a(i,mmaxp1) = a(i,mmaxp1) - 2.0_rk*alph*eb(ll,i)
                fk(i) = fk(i) - 2.0_rk*alph*eb(ll,i)*eb(ll,mmaxp1)
            end do
            a(mmaxp1,mmaxp1) = a(mmaxp1,mmaxp1) - 2.0_rk*alph
            fk(mmaxp1) = fk(mmaxp1) + 2.0_rk*alph*eb(ll,mmaxp1)
        end do
    end if

    do i = 1, mpnmax
        do k = i, mpnmax
            a(k,i) = a(i,k)
        end do
    end do

    do j = 1, Nmax
        fk(mmaxp1+j) =  - bb(j,mmaxp1)
    end do

    if (icops >= 2) fk(mpnmax) = value

    call gelg(fk, a, mpnmax, 1, 1e-7, ier)

    psicon = fk(mmaxp1)
    energy = 0.0_rk
    do i = 1, Mmax
        ki = i + 1
        if (ki <= Mmax) then
            do k = ki, Mmax
                energy = energy + cl(i,k)*fk(i)*fk(k)
            end do
        end if
    end do

    fk(mmaxp2) = energy

    return

end subroutine xcur