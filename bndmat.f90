!   -----------------------------------------------------------------------
!
!     Computes the Boundary Matrix for subsequent use in routine EQSIL
!     in rectangular torus cross-section with coordinates r,z:  
!
!     bndmat.f90
!-----------------------------------------------------------------------
      
subroutine bndmat()    
    use solve_mod

    real(rk), dimension(:), allocatable :: rt, zt

    real(rk) :: r0, ra, rb, za, zb, arh, azh, ar, az, gfl
    integer(ik) :: i, j, j1, j2, kp, m2, n2, nl


    allocate(rt(llp), zt(llp))

    ar = 2.0_rk/real(Mm1, rk)
    r0 = 1.0_rk/alpha
    az = sh*ar

! 
! Index vector interior neighbours to boundary, bottom-top
!
    j1 = Mr
    j2 = Mr*(Nz-2)
    kp = 1
    do i = 2, Mm1
        ip(kp) = i + j1
        ip(kp+1) = i + j2
        kp = kp + 2
    end do
!
! Left - Right
!
    do j = Mr, j2, Mr
        ip(kp) = j + 2
        ip(kp+1) = j + Mm1
        kp = kp + 2
    end do
    kp = kp - 1
!
!  r - coordinates of boundary points, indexvector
!
    ra = r0 - 1.0_rk      
    rb = r0 + 1.0_rk  
    za = 0.0_rk  
    zb = float(Nm1)*az
    m2 = Mr - 2
    j2 = Mr*Nm1

!
! Bottom - Top
!
    nl = 1
    do i = 1, m2
        rt(nl) = ra + float(i)*ar
        rt(nl + 1) = rt(nl)
        zt(nl) = za
        zt(nl + 1) = zb
        jp(nl) = i + 1
        jp(nl + 1) = i + 1 + j2
        nl = nl + 2      
    end do

!
! Left - Right
!
    n2 = Nz - 2
    do j = 1, n2
       rt(nl) = ra
       rt(nl + 1) = rb
       zt(nl) = float(j)*az
       zt(nl + 1) = zt(nl)
       jp(nl) = j*Mr + 1
       jp(nl + 1) = (j + 1)*Mr
       nl = nl + 2
    end do

!
! Matrix elements
!
    arh = 0.5_rk*ar
    azh = 0.5_rk*az
    do i = 1, kp 
!
! Bottom - Top
!
        nl = 1
        do j = 1, m2
            aux(nl,i) = gfl(rt(i), rt(nl), zt(i) - zt(nl), ar)/(sh*rt(nl)) 
            aux(nl + 1, i) = gfl(rt(i), rt(nl + 1), zt(i) - zt(nl + 1), ar)/(sh*rt(nl + 1))
            nl = nl + 2
        end do
!
! Left - Right
!
        do j = 1, n2
            aux(nl,i) = gfl(rt(i),rt(nl),zt(i) - zt(nl), az)*sh/(rt(nl) + arh)
            aux(nl + 1,i) = gfl(rt(i),rt(nl + 1),zt(i) - zt(nl + 1), az)*sh/(rt(nl + 1) - arh)
            nl = nl + 2   
        end do
    end do

    return        
        
end subroutine bndmat

!-----------------------------------------------------------------------

function gfl (rv, rst, zv, del)
    use solve_mod

    real(rk) :: ak, x, xdl, rv, rst, zv, del, gfl

    ak = 4.0_rk*rv*rst/((rv + rst)**2 + zv**2)
    x = ((rv - rst)**2 + zv**2)/((rv + rst)**2 + zv**2)   
    if (x.eq.0.) then  
        xdl = 2.0_rk*(log(del/(4.0_rk*rv)) - 1.0_rk)
     else 
        xdl = log(x)
     end if      

    gfl = sqrt(rv*rst/ak)*((1.0_rk - 0.5_rk*ak)*k(x) - e(x))/pi

    return
    
    contains
! 
!   Internal  functions
!

    function p1(x)
        real(rk) x, p1

        p1 = (((.01736506451_rk*x + .04757383546_rk)*x + .06260601220_rk)*x + .44325141463_rk)*x + 1.0_rk
        
    end function p1
  ! 
    function p2(x)
        real(rk) x, p2
    
        p2 = (((.00526449639_rk*x + .04069697526_rk)*x + .09200180037_rk )*x + .24998368310_rk)*x
    
    end function p2
  ! 
    function p3(x)
        real(rk) x, p3
        
        p3 = (((.01451196212_rk*x + .03742563713_rk)*x + .03590092383_rk)*x + .09666344259_rk)*x + 1.38629436112_rk
    
    end function p3
  ! 
    function p4(x)
        real(rk) x, p4
    
        p4 = (((.00441787012_rk*x + .03328355346_rk)*x + .06880248576_rk)*x + .12498593597_rk)*x + .5_rk
    
    end function p4
  ! 
    function e(x)
        real(rk) x, e
        
        e = p1(x) - p2(x)*xdl
        
    end function e
  
    function k(x)
        real(rk) x, k

        k = p3(x) - p4(x)*xdl
    end function k
  ! 

end function gfl   
 
