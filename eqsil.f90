subroutine eqsil(q)
    use solver_mod

   real(rk), dimension(*) :: q
   real(rk), dimension (:), allocatable :: qq
   integer(ik) :: i, k, l, npn
   real(rk) :: sum

   allocate (qq(MN))

   qq(1:MN) = q(1:MN)

   npn = Nm1*Mr
   do i = 1, Mr
      qq(i) = 0.0_rk
      qq(i + npn) = 0.0_rk
   end do
   do k = 1, MN, Mr
      qq(k) = 0.0_rk
      qq(k + Mm1) = 0.0_rk
   end do

   call flux(qq)

   do i = 1, llp
      sum = 0.0_rk
      do l = 1, llp
         sum = sum + qq(ip(l))*aux(l,i)
      end do
         q(jp(i)) = q(jp(i)) + sum
   end do

   call flux(q)

   deallocate (qq)
   
   return 

end subroutine
