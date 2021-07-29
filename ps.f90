function ps (rpa, zpa, rp, zp, l)

    use solver_mod

    real(rk), intent(in) :: rpa, zpa, rp, zp

    real(rk),dimension(3,2) :: sigc
    real(rk),dimension(2) :: sigp, c, d
    real(rk),dimension(3) :: a,b
    real(rk) acl,args,alg,abl,resk,rese, ps
    integer :: k,l

    data sigc /1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 0.0_rk/ sigp /1.0_rk, -1.0_rk/ 
    data a /1.3862944_rk, .1119723_rk, .0725296_rk/ b /.5_rk, .1213478_rk, .0288729_rk/ 
    data c /.4630151_rk, .1077812_rk/ d /.2452727_rk, .0412496_rk/
!
    ps = 0.0_rk
    do k = 1, 2
        acl = 1.0d-06
        args = 4.0_rk*rpa*rp/((rpa + rp)**2 + (zp-sigp(k)*zpa)**2)
        abl = 1.0_rk - args
        if (abl .ge. 1.0d-06) acl = abl
        alg = log(1.0_rk/acl)
        resk = (a(3)*acl + a(2))*acl + a(1) +  alg*((b(3)*acl + b(2))*acl + b(1))
        rese = (c(2)*acl + c(1))*acl + 1.0_rk + alg*(d(2)*acl + d(1))*acl
        ps = ps + sigc(l,k)*sqrt(rpa*rp/args)*((1.0_rk - 0.5_rk*args)*resk - rese)/pi
    end do

    print *, 'ps = ', ps

    return 

      
end function