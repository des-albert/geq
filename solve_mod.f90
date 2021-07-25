module solve_mod

    implicit none

    integer, parameter :: single = kind(1.0)
    integer, parameter :: double = kind(1.0d0)
    integer, parameter :: quad = selected_real_kind(p=30)
    integer, parameter :: rk = double

    integer, parameter :: int_def = kind(1)
    integer, parameter :: ik = int_def
    integer, parameter :: Mc = 15
    integer, parameter :: Ng = 1
    integer, parameter :: Nmax = 6
    integer, parameter :: llmax = 16

    integer(ik) :: Nexp, Mr, Nz, MN, Mm1, Nm1, llp, mpnmax
    integer(ik) :: jdes, jaxis, ndes, naxis, icycle, idecis, idol 
    integer(ik) :: Mmax, irsp, izsp
    integer, dimension(6) :: ityp

    integer, dimension(:), allocatable :: ip, jp

    real(rk), PARAMETER :: pi = 3.141592653589793238462643383279502884197_rk

    real(rk) :: Rmin, Rmax, Zmin, Zmax, Rmpl, Zmpl, Apl, error, alpha, sh, ss
    real(rk) :: dr, dz, alph, totcurr, betapol, alfac, psicon, fmaxa, fabs, value
    real(rk) :: c0pr, c0btr, xind, pint, pintvo, bt2int, betapvo, betat, betap, beshysu, beshyvo
    real(rk), dimension(:), allocatable :: R, Z, com, cjt, pr, bt2
    real(rk), dimension(:,:), allocatable :: aux
    real(rk), dimension(Nmax, Mc) :: bb
    real(rk), dimension(llmax, Mc) :: eb
    real(rk), dimension(Mc + 1, Mc) :: cl
    real(rk), dimension(Mc + Nmax) :: fk

!   Resize g and f when Nexp > 6

    real(rk), dimension(65, 65) :: g
    real(rk), dimension(4225) :: f

    equivalence (f(1), g(1,1))

end module





  