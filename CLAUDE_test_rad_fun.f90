! Small driver for the spheroidal_radial module.
! Build:  gfortran -O2 spheroidal_radial.f90 test_rad_fun.f90 -o test_rad_fun
program test_rad_fun
   use spheroidal_radial, only: dp, rad_fun
   implicit none
   integer, parameter :: ne = 5
   complex(dp) :: r1f(ne), r1d(ne), r2f(ne), r2d(ne)
   integer :: kob, m, n
   complex(dp) :: c
   real(dp) :: xi

   ! prolate, m = 1, c = 5, xi = 1.5
   kob = 0; m = 1; c = (5.0_dp, 0.0_dp); xi = 1.5_dp
   call rad_fun(kob, m, ne, c, xi, 1.0e-15_dp, r1f, r1d, r2f, r2d)
   call show('prolate')

   ! oblate, m = 1, c = 5 + 0.5i, xi = 0.5
   kob = 1; m = 1; c = (5.0_dp, 0.5_dp); xi = 0.5_dp
   call rad_fun(kob, m, ne, c, xi, 1.0e-15_dp, r1f, r1d, r2f, r2d)
   call show('oblate')

contains

   subroutine show(title)
      character(*), intent(in) :: title
      print '(/a,"  m =",i3,"  c =",2f8.3,"  xi =",f8.3)', title, m, c, xi
      print '(a4,4(2x,a27))', 'n', 'R1', 'dR1/dxi', 'R2', 'dR2/dxi'
      do n = 1, ne
         print '(i4,4(2x,"(",es12.5,",",es12.5,")"))', m + n - 1, &
            r1f(n), r1d(n), r2f(n), r2d(n)
      end do
   end subroutine show

end program test_rad_fun
