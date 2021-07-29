!-----------------------------------------------------------------------
!     solver.f90
!
!     Author : Des Albert
!
!     Lackner Tokamak Equilibrium Solver
!
!
!     Modified by Des Albert
!
!------------------------------------------------------------------------------- 
! 
!    Data input Lines 
!  
! 1  Title
! 2  rmin,rmax,zmin,zmax,err,meshfg,mprfg
!                 min and max r,z values of calculation grid.
!                 err is iteration tolerence, usually 0.0001
!                 mprfg = 0 no output of progressive iteration data
! 3  Blank
! 4  rmpl, offset, apl, el, tri, rxpsn, zxpsn
!
! 5  Blank
! 6  rac,zac,exc
!                 Groups of Poloidal Field Conductors.
!                 Each group can contain several line, separated by 0.
!                 Final group terminated by two 0. cards.
!                 rac,zac are r,z coordinates of pf coils.
!                 exc is current multiplier.
!                 Each group treated as one variable. Individual
!                 currents = overall x exc. rlc is not used.
!
! 7  icops,value  icops = 1  flux and net ampereturns free
!                         2  net ampereturns of system = value
!                         3  flux at plasma boundary = value
! 8,9  Blank
! 10,11  current,betpr,alfac
!                 Total plasma current, approximate poloidal beta,
!                 current function factor
!  raxis,zaxis,zdes,alp  r,z coordinates of Magnetic Axis.
!                 
!                 This must be a grid point. If not is moved to the
!                 nearest. These values alters as the program iterates.
!                 alp is the balance factor for the optimisation process
!                 between a good fit to the plasma points with ix = 4,
!                 and a low pf coil cost function. High alp implies good
!                 fit to points. alp usually in range 10**15-10**5.
!                 Present pf coil cost function sig(i**2).
! 
!  Note:
!  psi(Wb)  =  0.8*pi*pi*psi(code)
!  currents are in  MA
! 

program solver
    use solver_mod

    integer(ik) :: meshfg, mprfg, nof, jn, ll, kk
    integer(ik) :: i, j, k, ii, icm, lp1, icl, nlines, na
    integer(ik), dimension(Mc) :: ic
    integer(ik), dimension(30) :: np

    real(rk) :: Offset, El, Rxpsn, Zxpsn, rac, zac, exc, tri, trixp, elxp, ang, anga, angb, ang1, ang2
    real(rk) :: rc1, rc2, al1, al2, gfl, raxis, zaxis, zdes, alp, xt1, xt2, xt3, curr 
    real(rk) :: dzh, qsc, vdc, btvac, br1, fb, dov 
    real(rk), dimension(Ng, Mc) :: Ra, Za, Ex, Rl
    real(rk), dimension(16) :: Rcc, Zcc
    real(rk), dimension(Nmax) :: Rc, Zc
    real(rk), dimension(:), allocatable :: expsi, fool
    real(rk), dimension(:,:), allocatable :: psiext
    real(rk),dimension(30) :: psiq, rmn, rmx, qs, vd
    real(rk),dimension(100) :: zh, qdv 

    Nexp = 6    

    Mr = 1 + 2**Nexp   
    Nz = Mr
    MN = Mr * Nz
    Mm1 = Mr - 1
    Nm1 = Nz - 1
    llp = 2*(Mr + Nz) - 8

!   Read data file    

    open (8, file='solver.dat')
    read(8,*)
    read(8,*) Rmin, Rmax, Zmin, Zmax, error, meshfg, mprfg
    read(8,*)
    read(8,*) Rmpl, Offset, Apl, El, tri, Rxpsn, Zxpsn

    if (meshfg.gt.0) then
        Rmin = Rmpl - 32.0_rk*Apl/20.0_rk
        Rmax = Rmpl + 32.0_rk*Apl/20.0_rk
        Zmin = Offset - 32.0_rk*(Offset + abs(Zxpsn))/20.0_rk
        Zmax = Offset + 32.0_rk*(Apl*El)/20.0_rk
      end if


    allocate(R(Mr), Z(Nz), cjt(Mr), pr(Mr), bt2(Mr))
    allocate (ip(llp), jp(llp)) 
    allocate (aux(llp, llp))

    dr = (Rmax - Rmin)/real(Mm1, rk)
    dz = (Zmax - Zmin)/real(Nm1, rk)

    do i= 1, Mr
        R(i) = Rmin + (i - 1)*dr
    end do
    do j = 1, Nz
        Z(j) = Zmin + (j - 1)*dz
    end do

    alpha = (Rmax - Rmin)/(Rmax + Rmin)
    sh = dz/dr
    ss = sh**2

    write (*,'(a,f7.3,a,f7.3,a,f7.3)') "Rmin  = ", Rmin, " Rmax = ", Rmax," dr = ", dr 
    write (*,'(a,f7.3,a,f7.3,a,f7.3)') "Zmin  = ", Zmin, " Zmax = ", Zmax," dz = ", dz
    write (*,'(a,f7.3,a,f7.3,a,f7.3)') "alpha = ", alpha, "   sh = ", sh," ss = "
!    
! Boundary Matrix Setup
!
    call bndmat()

!
! Read in Poloidal Field Coil Data
!   

    k = 0
    read (8,*)
    do while (.true.)
       k = k + 1
       ic(k) = 0
       do while (.true.)
          read (8,*) rac, zac, exc
          if (rac .le. 0.0_rk) go to 10
          ic(k) = ic(k) + 1       
          Ra(ic(k),k) = rac
          Za(ic(k),k) = zac
          Ex(ic(k),k) = exc
          Rl(ic(k),k) = 1.0e-20*rac
       end do

 10    if (ic(k) .lt. 1) go to 20
    end do
20  Mmax = k - 1
    if (mprfg.ne.0) then
        write(*,'(a)') "Conductor groups available for optimization"
        do k = 1, Mmax
            write (*,'(a,i3)') "Group :", k
            do i = 1, ic(k)
                write (*,'(5x,f7.3,3x,f7.3,3x,f7.3)') Ra(i,k), Za(i,k), Ex(i,k)
            end do
        end do
    end if 

!
!  Conditions to be satisfied by resulting equilibrium
!
    elxp = abs(Offset - Zxpsn)/Apl
    trixp = (Rmpl - Rxpsn)/Apl
!
    do j = 1, 3
       ityp(j) = j
       Rc(j) = Rxpsn
       Zc(j) = Zxpsn
    end do

    ityp(4) = 1
    Rc(4) = Rmpl - Apl
    Zc(4) = Offset
    ityp(5) = 1
    Rc(5) = Rmpl + Apl
    Zc(5) = Offset
    ityp(6) = 1
    Rc(6) = Rmpl - Apl*tri
    Zc(6) = Offset + Apl*El
    do j = 1, 4
       ang = j*pi/10.0_rk
       Rcc(j) = Rmpl + Apl*cos(ang + tri*sin(ang))
       Zcc(j) = Offset + El*apl*sin(ang)
       Rcc(j+4) = Rmpl + Apl*cos(ang + pi/2.0_rk + tri*sin(ang + pi/2.0_rk))
       Zcc(j+4) = Offset + El*Apl*sin(ang + pi/2.0_rk)
    end do
    
    do j = 1, 4
        al1 = Apl*(((1.0_rk + trixp)*(1.0_rk + trixp)) + elxp*elxp)/(2.0_rk*(1.0_rk + trixp))
        al2 = Apl*(((1.0_rk - trixp)*(1.0_rk - trixp)) + elxp*elxp)/(2.0_rk*(1.0_rk - trixp))
        anga = atan(2.0_rk*elxp*(1.0_rk + trixp)/(elxp*elxp - (1.0_rk + trixp)*(1.0_rk + trixp)))
        angb = atan(2.0_rk*elxp*(1.0_rk - trixp)/(elxp*elxp - (1.0_rk - trixp)*(1.0_rk - trixp)))
        rc1 = Rmpl + Apl - al1
        rc2 = Rmpl - Apl + al2
        ang1 = anga*j/5.0_rk
        ang2 = angb*j/5.0_rk
        Rcc(j+8) = rc1 + al1*cos(ang1)
        Zcc(j+8) = Offset - al1*sin(ang1)
        Rcc(j+12) = rc2 - al2*cos(ang2)
        Zcc(j+12) = Offset - al2*sin(ang2)
    end do

    if (mprfg.ne.0) then
        write (*,*)
        write (*,*) '  Single null case: boundary points '

        do j = 1, 16
            write (*,'(5x,i4,2x,f8.3,2x,f8.3)') j, Rcc(j), Zcc(j)
        end do

    end if 

    allocate(expsi(MN), fool(MN), psiext(MN, Mc))

    do kk = 1, Mmax
        icl = ic(kk)
        do i = 1, Nz
           nof = (i - 1)*Mr
           do j = 1, Mr
              jn = nof + j
              expsi(jn) = 0.0_rk
           end do
        end do
       
        do i = 1, icl
            if (.not. (((Zmax - Za(i,kk))*(Zmin - Za(i,kk)).le.0.0_rk) .and.  &
                ((Rmax - Ra(i,kk))*(Rmin - Ra(i,kk)) .le. 0.0_rk))) then
              do k = 1, Nz, Nm1
                 nof = (k - 1)*Mr
                 do j = 1, Mr
                    jn = nof + j
                    expsi(jn) = expsi(jn) + Ex(i,kk)*gfl(R(j), Ra(i,kk), Z(k) - Za(i,kk), 0.)
                 end do
              end do

              do k = 1, Nz
                 nof = (k - 1)*Mr
                 do j = 1, Mr, Mm1
                    jn = nof + j                    
                    expsi(jn) = expsi(jn) + Ex(i,kk)*gfl(R(j), Ra(i,kk), Z(k) - Za(i,kk), 0.)
                 end do
              end do

           end if
        end do

        call eqsil(expsi)

        do i = 1, icl
            if ((.not. (((Zmax - Za(i,kk))*(Zmin - Za(i,kk)) .gt. 0.0_rk) .or. &
                    ((Rmax - Ra(i,kk))*(Rmin - Ra(i,kk)) .gt. 0.0_rk)))) then
                do k = 1, Nz
                   nof = (k - 1)*Mr
                   do j = 1, Mr
                      jn = nof + j                    
                      expsi(jn) = expsi(jn) + Ex(i,kk)*gfl(R(j), Ra(i,kk), Z(k) - Za(i,kk), 0.)
                   end do
                end do
             end if
          end do

        do j = 1, MN
            psiext(j,kk) = expsi(j)
        end do

!
!  Computation of matrix elements for exact conditions
!
        call splnco (expsi)
        do j = 1, Nmax
            call condit(expsi, Rc(j), Zc(j), ityp(j), bb(j,kk))
        end do
        do k = 1, llmax
            call condit(expsi, Rcc(k), Zcc(k), 1, eb(k,kk))
        end do

!
!  Computation of inductances
!    

        cl(kk,kk) = 0.0_rk
        cl(Mmax+1, kk) = 0.0_rk
        do i = 1, icl
            cl(Mmax + 1, kk) = cl(Mmax + 1,kk) + Ex(i,kk)
            cl(kk,kk) = cl(kk,kk) + Ex(i,kk)**2*1.0e6*(0.58_rk + log(Ra(i,kk)/rl(i,kk)))/(2.0_rk*pi)
        end do

        do i = 1, icl
            ii = i + 1
            if (ii .le. ic(kk)) then
                do j = ii, icl
                    cl(kk,kk) = cl(kk,kk) + 2.0_rk*Ex(i,kk)*Ex(j,kk)*gfl(Ra(j,kk), Ra(i,kk), Za(j,kk)- Za(i,kk),ar)
                end do
            endif
        end do
        
        lp1 = kk + 1
        if (lp1 .le. Mmax) then
            do k = lp1, Mmax
                icm = ic(k)
                cl(kk,k) = 0.0_rk
                do i = 1, icl
                    do j = 1, icm
                        cl(kk,k) = cl(kk,k) + Ex(i,kk)*Ex(j,k)*gfl(Ra(j,k), Ra(i,kk), Za(j,k) - Za(i,kk),ar)
                    end do
                end do  
            end do
        endif

    end do

!
! Computation of a new case
!

    read(8,*)
    read(8,*) icops, value    
    mpnmax = Mmax + Nmax + 1
    if (icops.ge.2) then
        mpnmax = mpnmax + 1
    end if
    read(8,*)
    read(8,*) totcurr, betapol, alfac 
    read(8,*)
    read(8,*) raxis, zaxis, zdes, alp
    jdes = int(0.1 + (raxis - R(1))/dr) + 1
    jaxis = int(0.1 + (raxis - R(1))/dr) + 1
    raxis = R(jaxis)
    naxis = int(0.1 + (zaxis - Z(1))/dz) + 1
    zaxis = Z(naxis)
    ndes = int(0.1 + (abs(zdes) - Z(1))/dz) + 1
    if (zdes .gt. 0.0_rk) then
        zdes = z(ndes)
    end if
    write(*,'(a,f8.3,a,f8.3)') 'Magnetic Axis  r = ',raxis, ' z = ', zaxis
    write (*,'(a,f8.3,a,e10.3)') 'Rail limiter   z = ', zdes, ' alp factor = ',alp
    if (llmax.gt.0) then
        alph = alp*2.0_rk*pi/(llmax*raxis)
    end if

    close(8)

    call startt()

!
! Begin Iterations
!
    icycle = 1
    idecis = 0

    allocate(com(Mr))

    do while (icycle.le.20)
        write(6,'(a,i2,a)') ' ==== Cycle number ', icycle,' ===='

        call eqsil(f)

        call compar()

        do j = 1, MN
            fool(j) = 0.0_rk
            expsi(j) = f(j)
        end do  
        
        call splnco(expsi)

        do j = 1, Nmax
            call condit(expsi, Rc(j), Zc(j), ityp(j), bb(j,Mmax+1))
        end do
        do ll = 1, llmax
            call condit(expsi, Rcc(ll), Zcc(ll), 1, eb(ll,Mmax+1))
        end do

        call xcur (expsi)
    
        xt1 = 0.0_rk
        xt2 = 0.0_rk
        xt3 = 0.0_rk
        do i = 1, Mmax
            do j = 1, ic(i)
                curr = Ex(j,i)*fk(i)
                xt1 = xt1 + curr**2
                xt2 = xt2 + abs(curr)
                xt3 = xt3 + abs(curr*Ra(j,i))
            end do
        end do
        if (mprfg.ne.0) write (*, '(a,f12.4,a,f12.4,a,f12.4)') ' SIG(I**2) = ',xt1,'  SIG(ABS(I)) = ', xt2, '  SIG(ABS(RI)) = ', xt3

        do k = 1, Mmax
            do j = 1,  MN
                    fool(j) = fool(j) + fk(k)*psiext(j,k)
                    f(j) = f(j) + fk(k)*psiext(j,k)
            end do
        end do

        call saddle()

        if (irsp .gt. 2) then
            if (mprfg.ne.0) write (6,'(a,f12.5,a,f12.5)') ' Saddle point r = ', R(irsp), ' z = ', Z(izsp)
        end if
        if (idecis.gt.0) GO TO 30

        call curf()

        icycle = icycle + 1 

    end do

30  call plotit()  

 

!
!       Magnetic and Geometric properties
!
    nlines = 20
    na = 20 
    na = na + nlines + 1
    psiq(1) = f(naxis*Mr + jaxis)
    fmaxa = psiq(1) - psicon
    zh(1) = psicon + fmaxa
    dzh = fmaxa/nlines
    DO l = 2, nlines + 1
        zh(l) = zh(l-1) - dzh
    end do
    zh(nlines + 1) = psicon
    DO i = nlines + 2, na
        zh(i) = zh(i-1) - dzh
    end do

 

    nqmax = 1
    rmn(1) = raxis
    rmx(1) = raxis
    do i = 1, na
        call contour (zh(i))
    end do

    goto 40 

    betap = betap*circf**2
    betapvo = betapvo*circf**2
    beshysu = betap/(0.4_rk*pi)
    beshyvo = betapvo/(0.4_rk*pi)

    vd(1) = vd(2)
    qs(1) = qs(2)
    np(1) = 0
    qsc = qs(1)
    vdc = vd(1)

    do l = 1,nqmax
        qdv(l) = qs(l)/vd(l)
    end do
    call smootl (qs, 7, br1, cr1)
    qs(1) = cr1
    call smootl (vd, 7, br1, cr1)
    vd(1) = cr1
    call smootl (qdv, 7, br1, cr1)
    btvac = sqrt(1.0_rk/qs(1)**2)
    dov = br1/(zh(2)-zh(1))
    do j = 1, nqmax
        fb = sqrt(btvac**2)
        qs(j) = fb*qs(j)
    end do
    bt2int = btvac**2*betat + bt2int
    betat = 2.*pint/bt2int

    write(*,'(5x,a,f12.5)') 'Area Integral of Pressure         ', pint
    write(*,'(5x,a,f12.5)') 'Volume Integral of Pressure       ', pintvo
    write(*,'(5x,a,f12.5)') 'Area Integral of B-tor**2         ', bt2int
    write(*,'(5x,a,f12.5)') 'Poloidal Beta (Lack.-Sur.)        ', betap
    write(*,'(5x,a,f12.5)') 'Poloidal Beta (Lack.-Vol.)        ', betapvo
    write(*,'(5x,a,f12.5)') 'Poloidal Beta (Shynia.-Sur.)      ', betashysu
    write(*,'(5x,a,f12.5)') 'Poloidal Beta (Shynia.-Vol.)      ', betashyvo
    write(*,'(5x,a,f12.5)') 'Toroidal Beta                     ', betat
    write(*,'(5x,a,f12.5)') 'Toroidal Vacuum Field             ', btvac
    write(*,'(5x,a,f12.5)') 'Plasma Internal Inductance        ', xind

40  continue

end program solver