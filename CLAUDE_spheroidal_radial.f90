!=============================================================================
!  spheroidal_radial  --  radial spheroidal functions R^(1)_mn, R^(2)_mn
!
!  Modern Fortran (2008+) translation of rad_fun.for (Fortran 77, Nov 2005).
!  The numerical algorithms are kept exactly as in the original; the changes
!  are structural (module instead of COMMON blocks, IMPLICIT NONE, INTENTs,
!  no arithmetic IF / computed GOTO / PAUSE / Hollerith, shared code factored
!  into helper functions).  See the notes at the end of this file.
!
!  Public entry point:
!     call rad_fun(kob, m, ne, c2, ksi0, eps, r1f, r1d, r2f, r2d)
!
!     kob  - 0: prolate, 1: oblate functions
!     m    - index m
!     ne   - maximum index n  (functions for n = m, ..., m+ne-1)
!     c2   - argument c (complex)
!     ksi0 - argument xi (real)
!     eps  - accepted for compatibility; NOT used (the original code
!            overwrote it with 1e-15, which is what is used here as well)
!
!  Results (complex(dp), size ne):
!     r1f - R^(1)_mn(c,xi),  r1d - dR^(1)_mn/dxi
!     r2f - R^(2)_mn(c,xi),  r2d - dR^(2)_mn/dxi
!=============================================================================
module spheroidal_radial
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: dp, rad_fun

   integer, parameter :: nterms = 330
   real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
   complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)   ! former S
   complex(dp), parameter :: ci = (0.0_dp, 1.0_dp)      ! former S1

   ! ---- former COMMON /K1/ -------------------------------------------------
   real(dp) :: aksi_g = 0.0_dp   ! xi^2 -/+ 1
   real(dp) :: ak = 1.0_dp       ! +1 prolate, -1 oblate
   integer  :: kob = 0           ! 0 prolate, 1 oblate (former K)
   integer  :: nk = 0
   ! ---- former COMMON /EPS3/, /EPS33/, /FACT/, /LEG/ -------------------------
   real(dp) :: eps3 = 1.0e-200_dp
   real(dp) :: eps33 = 1.0e-200_dp
   real(dp) :: fact(300) = 0.0_dp          ! fact(k+1) = k!  (scaled for k >= 169)
   real(dp), dimension(2*nterms) :: p_leg = 0.0_dp, pd_leg = 0.0_dp, &
                                    q_leg = 0.0_dp, qd_leg = 0.0_dp
   ! ---- former COMMON /F2/ (parameters of the ODE right-hand sides) ----------
   complex(dp) :: f2_rl = czero, f2_c = czero
   integer     :: f2_m = 0, f2_i12 = 0
   ! ---- log file (former unit 7, 'sph_fun.out') ------------------------------
   integer :: lu = 0

   abstract interface
      subroutine rhs_fn(t, y, yp)
         import :: dp
         real(dp), intent(in) :: t
         complex(dp), intent(in) :: y(*)
         complex(dp), intent(out) :: yp(*)
      end subroutine rhs_fn
   end interface

contains

!=============================================================================
!  Radial spheroidal functions (main driver)
!=============================================================================
   subroutine rad_fun(kob_in, m, ne, c2, ksi0, eps, r1f, r1d, r2f, r2d)
      integer, intent(in) :: kob_in, m, ne
      complex(dp), intent(in) :: c2
      real(dp), intent(in) :: ksi0
      real(dp), intent(in) :: eps          ! not used, see header
      complex(dp), intent(out) :: r1f(ne), r1d(ne), r2f(ne), r2d(ne)

      character(*), parameter :: fmt_w = &
         '(1x,"ifun=",i3,2x,"L=",i3,2x,"W=",1pd13.6,2x,"W1=",d13.6,2x,"W2=",d13.6)'
      character(*), parameter :: fmt_ne = '(1x,"homFUNq NN>nterms",5x,"NN=",i5)'

      complex(dp) :: rlc2(nterms), rdc2(4*nterms+1), bdc2(4*nterms+1)
      complex(dp) :: rc1, rl, rvn, r1, r2, r3, r4, rr3, rr4
      real(dp) :: x, w, w0, w1, w2a, w2b, eps1, eps_lam
      integer :: ncc, ie, inum1, ifun1, l, lm, ifu, istep
      integer :: ier0, ier1, ier2, ier3, ier4

      call open_log()

      nk = ne + 40
      if (nk < 60) nk = 60
      if (ksi0 > 1.5_dp) nk = nk + 40
      if (nk > 4*nterms) nk = 4*nterms

      kob = kob_in
      if (kob /= 0) then
         ak = -1.0_dp
      else
         ak = 1.0_dp
      end if
      x = ksi0
      if (kob == 0) aksi_g = x**2 - 1.0_dp
      if (kob == 1) aksi_g = x**2 + 1.0_dp

      ! (the original stopped here with PAUSE; now only a warning is issued)
      if (kob == 0 .and. x < 1.00001_dp) then
         write (*, '(1x,3("!"),"    Prolate x > 1  x = ",f5.3)') x
         write (lu, '(1x,3("!"),"    Prolate x > 1  x = ",f5.3)') x
      end if
      if (kob == 1 .and. x < 0.00001_dp) then
         write (*, '(1x,3("!"),"    Oblate x > 0   x = ",f5.3)') x
         write (lu, '(1x,3("!"),"    Oblate x > 0   x = ",f5.3)') x
      end if

      eps_lam = 1.0e-15_dp
      eps1 = 1.0e-12_dp
      eps3 = 1.0e-200_dp
      eps33 = eps3
      call init_factorials()

      ! Wronskian reference value (real part, as in the original)
      w1 = real(1.0_dp/(c2*(ksi0**2 - 1.0_dp + 2*kob)), dp)
      ncc = int(real(c2, dp))
      if (ne > nterms) then
         write (lu, fmt_ne) ne
         write (*, fmt_ne) ne
         return
      end if
      rc1 = c2

      ! eigenvalues lambda
      call lambda(kob, m, ne, c2, eps_lam, rlc2, ie)
      if (ie /= 0) return

      inum1 = nk
      ifun1 = 11
      if (kob == 0) then
         ifun1 = 44
         if (ksi0 >= 1.5_dp) ifun1 = 22
      end if
      if (kob == 1) then
         ifun1 = 33
         if (ksi0 >= 1.5_dp) ifun1 = 22
      end if

      ! Legendre functions
      call funlegnn(m, ksi0, inum1)

      w0 = 0.0_dp
      r1 = czero; r2 = czero; r3 = czero; r4 = czero
      rr3 = czero; rr4 = czero

      do l = 1, ne
         lm = l + m - 1
         rl = rlc2(l)

         call cdcof4a(rl, rdc2, nk, bdc2, nk, 2, rvn, 1, m, lm, rc1, ier1)
         if (abs(bdc2(m+1)) <= eps33*1.0e20_dp) eps33 = eps33/1.0e20_dp

         if (kob == 1 .and. l > ncc .and. ksi0 > 0.25_dp .and. ksi0 < 0.8_dp) ifun1 = 11

         ! state machine replacing the original GOTO network (labels 11,13,22,33,44,5)
         istep = ifun1
         do
            select case (istep)

            case (11)   ! LEG/LEG
               call cdrf12cc(r1, r2, r3, r4, 1, m, lm, rc1, rdc2, nk, bdc2, nk, ier2)
               w = real(r1*r4 - r2*r3, dp)
               w2a = abs(w/w1 - 1.0_dp)
               if (w2a > eps1) write (*, fmt_w) ifun1, l, w, w1, w2a
               if (w2a > w0 .and. kob == 0) then
                  istep = 13
               else
                  istep = 5
               end if

            case (13)   ! LEG/DIFF (JAF/DIFF)
               if (kob == 0) call cdrsf20(rr3, rr4, ksi0, r1, r2, m, rc1, rl, ier4)
               w = real(r1*rr4 - r2*rr3, dp)
               w2b = abs(w/w1 - 1.0_dp)
               ifu = 13
               if (ifun1 == 22) ifu = 23
               if (ifun1 == 44) ifu = 43
               if (w2b > eps1) write (*, fmt_w) ifu, l, w, w1, w2b
               if (w2b < w0) then
                  r3 = rr3
                  r4 = rr4
               end if
               istep = 5

            case (22)   ! BES/BES
               call cdrb12cc(r1, r2, r3, r4, 1, m, lm, rc1, ksi0, rdc2, nk, ier3)
               w = real(r1*r4 - r2*r3, dp)
               w2a = abs(w/w1 - 1.0_dp)
               if (w2a <= eps1) then
                  istep = 5
               else if (kob == 0) then
                  istep = 13
               else
                  istep = 5
               end if

            case (33)   ! DIFF/DIFF
               call cdrsf212(r1, r2, r3, r4, 2, m, lm, rc1, ksi0, rl, &
                             rdc2, nk, bdc2, m+1, ier4)
               w = real(r1*r4 - r2*r3, dp)
               w2a = abs(w/w1 - 1.0_dp)
               if (w2a > eps1) then
                  write (*, fmt_w) ifun1, l, w, w1, w2a
                  istep = 11
               else
                  istep = 5
               end if

            case (44)   ! JAF/JAF
               call cdrg1cn(r1, r2, r3, r4, 1, m, lm, rc1, ksi0, &
                            rdc2, nk, bdc2, nk, rl, ier0)
               w = real(r1*r4 - r2*r3, dp)
               w0 = abs(w/w1 - 1.0_dp)
               if (w0 > eps1) then
                  write (*, fmt_w) ifun1, l, w, w1, w0
                  istep = 13
               else
                  istep = 5
               end if

            case default   ! label 5: final results
               exit
            end select
         end do

         r1f(l) = r1
         r1d(l) = r2
         r2f(l) = r3
         r2d(l) = r4

         write (lu, *) 'n ', l
         write (lu, *) 'R^(1) ', r1
         write (lu, *) 'R^(1)" ', r2
         write (lu, *) 'R^(2) ', r3
         write (lu, *) 'R^(2)" ', r4
         w = real(r1*r4 - r2*r3, dp)
         w2a = abs(w/w1 - 1.0_dp)
         write (lu, *) 'W, log W ', w2a, log10(w2a + 1.0e-100_dp)
         write (lu, '(1x,61("."))')
      end do
   end subroutine rad_fun

!-----------------------------------------------------------------------------
   subroutine open_log()
      logical :: opened
      inquire (file='sph_fun.out', opened=opened, number=lu)
      if (.not. opened) then
         open (newunit=lu, file='sph_fun.out', status='unknown', &
               position='append', action='write')
      end if
   end subroutine open_log

!-----------------------------------------------------------------------------
   subroutine init_factorials()
      integer :: i, j
      fact(1) = 1.0_dp
      fact(2) = 1.0_dp
      do i = 3, 170
         fact(i) = fact(i-1)*(i - 1.0_dp)
      end do
      fact(170) = fact(170)*1.0e-300_dp
      do j = 170, 297
         fact(j+1) = j*fact(j)
      end do
   end subroutine init_factorials

!=============================================================================
!  Lambda -- eigenvalues for complex c (continuation from real c if needed)
!
!  The original contained a cache of eigenvalues for previous calls that was
!  switched off (icc = 1 hard-coded); that unreachable code has been removed.
!=============================================================================
   subroutine lambda(k, m, ne, c2, eps, rlc2, ie)
      integer, intent(in) :: k, m, ne
      complex(dp), intent(in) :: c2
      real(dp), intent(in) :: eps
      complex(dp), intent(inout) :: rlc2(nterms)
      integer, intent(out) :: ie

      complex(dp) :: cc2a, rlcc2(nterms)
      real(dp) :: stepr, stepi, stepi0, cc2, imc2, ac2, aimc2, epss

      ie = 0
      cc2 = real(c2, dp)
      imc2 = aimag(c2)

      ! (several values were assigned to stepi0 in the original and then
      !  unconditionally overwritten by this one)
      stepi0 = 0.0021_dp

      if (stepi0 > imc2) then
         stepi0 = imc2/3.0_dp
      else if (k == 0 .and. imc2 < 0.5_dp) then
         call cdlamp(rlc2, m, ne, c2, eps, ie)
         if (ie == 0) return
         print *, c2, ie
      else if (k == 1) then
         ie = 10
         if (imc2 < 0.11_dp .and. cc2 <= 22.0_dp) call cdlamo_s(rlc2, m, ne, c2, eps, ie)
         if (imc2 < 0.3_dp .and. cc2 > 22.0_dp) call cdlamo_l(rlc2, m, ne, c2, eps, ie)
         if (ie == 0) return
         if (ie /= 10) print *, c2, ie
      end if

      ! ---- label 484: continuation in Im(c) starting near the real axis -----
      stepr = 0.0_dp
      stepi = stepi0
      ac2 = 0.091_dp
      if (imc2 < 0.1_dp) ac2 = 0.0_dp
      if (k == 1) ac2 = 0.0_dp
      cc2a = cmplx(cc2, ac2, dp)

      if (k == 0) then
         do
            call cdlamp(rlcc2, m, ne, cc2a, eps, ie)
            if (ie == 0) exit
            cc2a = cmplx(cc2, 0.0_dp, dp)
         end do
      end if

      if (k == 1) then
         if (real(cc2a, dp) <= 22.0_dp) call cdlamo_s(rlcc2, m, ne, cc2a, eps, ie)
         if (real(cc2a, dp) > 22.0_dp) call cdlamo_l(rlcc2, m, ne, cc2a, eps, ie)
         if (abs(c2) - abs(cc2a) < 1.0e-5_dp) then
            rlc2(1:ne) = rlcc2(1:ne)
            return
         end if
      end if

      do   ! label 1444
         cc2a = cc2a + cmplx(stepr, stepi, dp)
         aimc2 = imc2 - aimag(cc2a)
         if (abs(real(cc2a, dp) - cc2) < stepr .or. aimc2 < 0.0_dp) cc2a = c2

         epss = eps
         if (k == 0) epss = 1.0e-10_dp
         call cdlamn(rlcc2, rlc2, m, ne, cc2a, epss, ie)
         if (ie /= 0) print *, cc2a, ie

         if (ie /= 0) then
            cc2a = cc2a - 2.0_dp*cmplx(stepr, stepi, dp)
            if (aimag(cc2a) < 0.0_dp) cc2a = cmplx(real(cc2a, dp), 0.0_dp, dp)
            stepi = stepi/2.1_dp
            if (stepi < 1.0e-6_dp) return
            cycle
         end if

         rlcc2(1:ne) = rlc2(1:ne)
         if (abs(c2) - abs(cc2a) < 1.0e-5_dp) exit
      end do
   end subroutine lambda

!=============================================================================
!  BESSJJ0 -- spherical Bessel functions j_n by backward recurrence
!  (num may be reduced on overflow, as in the original)
!=============================================================================
   subroutine bessjj0(a, num, besj)
      complex(dp), intent(in) :: a
      integer, intent(inout) :: num
      complex(dp), intent(inout) :: besj(*)

      complex(dp) :: b, c
      integer :: n, num1, i, i1
      logical :: restart

      do
         besj(num+1) = czero
         besj(num) = cmplx(1.0e-300_dp, 1.0e-300_dp, dp)
         n = 2*num + 1
         num1 = num - 1
         restart = .false.
         do i = 1, num1
            n = n - 2
            i1 = num - i
            besj(i1) = n*a*besj(i1+1) - besj(i1+2)
            if (abs(besj(i1)) > 1.0e300_dp) then
               num = i - 1
               restart = .true.
               exit
            end if
         end do
         if (.not. restart) exit
      end do

      n = 2*(num/2)
      b = 1.2533141373155002_dp*sqrt(a)
      c = b*besj(1)
      do i = 3, n, 2
         b = b*(i - 0.5_dp)*(i - 2.0_dp)/(i - 2.5_dp)/(i - 1.0_dp)
         c = c + b*besj(i)
      end do
      c = 1.0_dp/c
      do i = 1, num
         besj(i) = c*besj(i)
      end do
   end subroutine bessjj0

!=============================================================================
!  CESSEL0 -- spherical Bessel functions j_n and y_n (y_n with scaling)
!=============================================================================
   subroutine cessel0(a, num, besj, besy, iiy)
      complex(dp), intent(in) :: a
      integer, intent(inout) :: num
      complex(dp), intent(inout) :: besj(*), besy(*)
      integer, intent(inout) :: iiy(*)

      complex(dp), parameter :: factor = (1.0e250_dp, 0.0_dp)
      complex(dp) :: b, b1, b2, c
      integer :: n, num1, i, i1, i2
      logical :: restart

      do
         iiy(1:num+1) = 0
         besj(num+1) = czero
         besj(num) = cmplx(1.0e-300_dp, 1.0e-300_dp, dp)
         n = 2*num + 1
         num1 = num - 1
         restart = .false.
         do i = 1, num1
            n = n - 2
            i1 = num - i
            besj(i1) = n*a*besj(i1+1) - besj(i1+2)
            if (abs(besj(i1)) > 1.0e300_dp) then
               num = i - 1
               restart = .true.
               exit
            end if
         end do
         if (.not. restart) exit
      end do

      b1 = a*besj(1) - besj(2)
      b2 = -a*b1 - besj(1)
      n = 2*(num/2)
      b = 1.2533141373155002_dp*sqrt(a)
      c = b*besj(1)
      do i = 3, n, 2
         b = b*(i - 0.5_dp)*(i - 2.0_dp)/(i - 2.5_dp)/(i - 1.0_dp)
         c = c + b*besj(i)
      end do
      c = 1.0_dp/c
      do i = 1, num
         besj(i) = c*besj(i)
      end do

      ! Y
      besy(1) = -c*b1
      besy(2) = c*b2
      do i = 3, num
         i2 = i - 2
         if (iiy(i-1) == iiy(i2)) then
            besy(i) = (2.0_dp*i2 + 1.0_dp)*a*besy(i-1) - besy(i2)
         else
            besy(i) = (2.0_dp*i2 + 1.0_dp)*a*besy(i-1) - besy(i2)/factor
         end if
         iiy(i) = iiy(i-1)
         if (abs(besy(i)) > 1.0e300_dp) then
            besy(i) = besy(i)/factor
            iiy(i) = iiy(i-1) + 1
         end if
      end do
   end subroutine cessel0

!=============================================================================
!  cdcof4a -- expansion coefficients d_r (ad) and, for id = 2, the
!             coefficients for r < 0 (bd); for in /= 0 the normalisation vn
!=============================================================================
   subroutine cdcof4a(vl, ad, ia, bd, ib, id, vn, in, m, n, c, ier)
      integer, intent(in) :: ia, ib, id, in, m, n
      complex(dp), intent(in) :: vl, c
      complex(dp), intent(inout) :: ad(ia), bd(ib)
      complex(dp), intent(inout) :: vn
      integer, intent(out) :: ier

      complex(dp) :: rl(4*nterms), s2, rc2, rc4, ru, rb, rv, rw, rr
      complex(dp) :: snm, rf1, rf2, sg, rvn, adl
      real(dp) :: ac, eps4, am, amm, bm, f2, d1
      integer :: n1, i11, l1, l11, j, jj1, jj2, ii, ii1, ia1, l4, i, ij
      integer :: iff, il1, j1, j2, j3

      ac = real(c, dp)
      ier = 0
      eps4 = 1.0e-20_dp
      ! check of input parameters
      if (m < 0 .or. n < m) ier = 1
      if (ac < 0.0_dp) ier = 3
      if (id == 1 .and. ia <= ((n - m)/2 + 2)) ier = 4
      if (id == 0 .and. in /= 0) ier = 6
      if (ier /= 0) return
      am = m
      amm = 2.0_dp*m
      s2 = (1.0_dp, 0.0_dp)
      bm = 4.0_dp*am**2 - 1.0_dp
      rc2 = c*c
      rc4 = rc2*rc2
      n1 = n - m
      if (id == 0) return

      ! ---- expansion coefficients -------------------------------------------
      if (mod(n1+1, 2) == 0) then
         i11 = 1
      else
         i11 = 0
      end if
      l1 = (n1 - i11)/2
      l11 = l1 + 1
      rl(l11) = s2

      ! coefficients for r < n-m
      do j = 1, l1
         jj1 = n1 - 2*j
         ru = s2
         rb = (am + jj1)*(am + jj1 + 1.0_dp) + ak*rc2/2.0_dp*(1.0_dp - bm/(amm + 2.0_dp*jj1 - 1.0_dp)/ &
              (amm + 2.0_dp*jj1 + 3.0_dp)) - vl
         rv = s2/rb
         rw = rv
         jj2 = jj1
         if (jj2 > 1) then
            do
               ii = 2*m + jj2
               ii1 = ii + jj2
               rr = -jj2*(jj2 - 1.0_dp)*ii*(ii - 1.0_dp)*rc4/((ii1 - 1.0_dp)**2*(ii1 - 3.0_dp)* &
                    (ii1 + 1.0_dp))/rb
               rb = (am + jj2 - 2.0_dp)*(am + jj2 - 1.0_dp) + ak*rc2/2.0_dp*(1.0_dp - bm/(ii1 - 5.0_dp)/ &
                    (ii1 - 1.0_dp)) - vl
               rr = rr/rb
               ru = 1.0_dp/(1.0_dp + rr*ru)
               rv = rv*(ru - 1.0_dp)
               rw = rw + rv
               if (abs(rv/rw) <= eps4 .and. j >= l1/2) exit
               jj2 = jj2 - 2
               if (jj2 <= 1) exit
            end do
         end if
         rl(l11-j) = -(amm + jj1 + 1.0_dp)*(amm + jj1 + 2.0_dp)*ak*rc2/(amm + 2.0_dp*jj1 + 3.0_dp)/ &
                     (amm + 2.0_dp*jj1 + 5.0_dp)*rw
         rl(l11-j) = rl(l11-j)*rl(l11+1-j)
      end do

      ! coefficients for r > n-m
      ia1 = ia - 1
      do j = l11, ia1
         jj1 = n1 + 2*(j - l11) + 2
         ii = 2*m + jj1
         ii1 = ii + jj1
         ru = s2
         rb = (am + jj1)*(am + jj1 + 1.0_dp) + ak*rc2/2.0_dp*(1.0_dp - bm/(ii1 - 1.0_dp)/(ii1 + 3.0_dp)) &
              - vl
         rv = -jj1*(jj1 - 1.0_dp)*ii*(ii - 1.0_dp)*rc4/((ii1 - 1.0_dp)**2*(ii1 - 3.0_dp)* &
              (ii1 + 1.0_dp))/rb
         rw = rv
         jj2 = jj1
         do
            jj2 = jj2 + 2
            ii = 2*m + jj2
            ii1 = ii + jj2
            rr = -jj2*(jj2 - 1.0_dp)*ii*(ii - 1.0_dp)*rc4/((ii1 - 1.0_dp)**2*(ii1 - 3.0_dp)* &
                 (ii1 + 1.0_dp))/rb
            rb = (am + jj2)*(am + jj2 + 1.0_dp) + ak*rc2/2.0_dp*(1.0_dp - bm/(ii1 - 1.0_dp)/(ii1 + 3.0_dp)) &
                 - vl
            rr = rr/rb
            ru = 1.0_dp/(1.0_dp + rr*ru)
            rv = rv*(ru - 1.0_dp)
            rw = rw + rv
            if (abs(rv/rw) <= eps4) exit
         end do
         rl(j+1) = rw*(amm + 2.0_dp*jj1 - 1.0_dp)*(amm + 2.0_dp*jj1 + 1.0_dp)/((amm + jj1 - 1.0_dp)* &
                   (amm + jj1)*ak*rc2)
         rl(j+1) = rl(j)*rl(j+1)
         if (j > l11) rl(j+1) = rl(j+1)/rl(l11)
         if (abs(rl(j+1)) < eps3) then
            do i = j + 2, ia
               rl(i) = czero
            end do
            exit
         end if
      end do

      ! ---- coefficient d(n-m) -----------------------------------------------
      l4 = (n + m + i11)/2
      f2 = 1.0_dp/fact(l1+1)
      iff = 0
      do i = 1, l4
         f2 = f2*(l4 + i)
         if (abs(f2) > 1.0e300_dp) then
            f2 = f2*1.0e-300_dp
            iff = 1
            ! NOTE: the original used the (undefined) variable J3 here instead
            ! of L4; L4 is clearly what was intended.
            do ij = i, l4
               f2 = f2*(l4 + ij)
            end do
            exit
         end if
      end do
      il1 = 1
      if (mod(l1, 2) /= 0) il1 = -1
      d1 = il1*f2/2.0_dp**n1
      if (iff == 1) d1 = d1*1.0e300_dp
      snm = czero

      do j = 1, ia
         j1 = 2*(j - 1) + i11
         j2 = (j1 - i11)/2
         j3 = (j1 + 2*m + i11)/2
         rf2 = rl(j)/2.0_dp**j1
         iff = 0
         do i = 1, j3
            rf2 = rf2*(j3 + i)
            if (abs(rf2) > 1.0e300_dp) then
               rf2 = rf2*1.0e-300_dp
               iff = 1
               do ij = i, j3
                  rf2 = rf2*(j3 + ij)
               end do
               exit
            end if
         end do
         rf1 = rf2
         do i = 1, j2
            rf1 = rf1/i
         end do
         sg = (-1.0_dp)**j2*rf1
         if (iff == 1) sg = sg*1.0e300_dp
         if (abs(sg) < eps3 .and. j > ia/2) exit
         snm = snm + sg
      end do
      ad(l11) = d1/snm

      adl = ad(l11)
      do i = 1, ia
         ad(i) = rl(i)*adl
      end do

      ! ---- coefficients for r < 0 -------------------------------------------
      if (id == 2) then
         bd(1) = ad(1)
         do j = 2, ib
            jj1 = i11 - 2*(j - 1)
            ru = s2
            rb = (am + jj1)*(am + jj1 + 1.0_dp) + ak*rc2/2.0_dp*(1.0_dp - bm/(amm + 2.0_dp*jj1 - 1.0_dp)/ &
                 (amm + 2.0_dp*jj1 + 3.0_dp)) - vl
            rv = 1.0_dp/rb
            rw = rv
            jj2 = jj1
            do
               ii = 2*m + jj2
               ii1 = ii + jj2
               rr = -jj2*(jj2 - 1.0_dp)*ii*(ii - 1.0_dp)*rc4/((ii1 - 1.0_dp)**2*(ii1 - 3.0_dp)* &
                    (ii1 + 1.0_dp))/rb
               rb = (am + jj2 - 2.0_dp)*(am + jj2 - 1.0_dp) + ak*rc2/2.0_dp*(1.0_dp - bm/(ii1 - 5.0_dp)/ &
                    (ii1 - 1.0_dp)) - vl
               rr = rr/rb
               ru = 1.0_dp/(1.0_dp + rr*ru)
               rv = rv*(ru - 1.0_dp)
               rw = rw + rv
               if (abs(rv/rw) <= eps4) exit
               jj2 = jj2 - 2
            end do
            if (jj1 == -2*m - 2 .or. jj1 == -2*m - 1) then
               bd(j) = ak*rc2/(amm + 2.0_dp*jj1 + 3.0_dp)/(amm + 2.0_dp*jj1 + 5.0_dp)*rw
               if (i11 == 1) bd(j) = -bd(j)
            else
               bd(j) = -(amm + jj1 + 1.0_dp)*(amm + jj1 + 2.0_dp)*ak*rc2/(amm + 2.0_dp*jj1 + 3.0_dp)/ &
                       (amm + 2.0_dp*jj1 + 5.0_dp)*rw
            end if
            bd(j) = bd(j)*bd(j-1)
            if (abs(bd(j)) < eps33) then
               do i = j + 1, ib
                  bd(i) = czero
               end do
               exit
            end if
         end do
      end if

      ! ---- normalisation factor ---------------------------------------------
      if (in == 0) return
      vn = czero
      iff = 0
      if (abs(ad(1)) > 1.0e150_dp) iff = 1
      do j = 1, ia
         jj1 = 2*(j - 1) + i11
         rf1 = (1.0_dp, 0.0_dp)
         if (iff == 1) rf1 = (1.0e-300_dp, 0.0_dp)
         rf1 = rf1*ad(j)/(amm + 2.0_dp*jj1 + 1.0_dp)*ad(j)
         if (m /= 0) then
            do i = 1, 2*m
               rf1 = rf1*(jj1 + i)
            end do
         end if
         rvn = rf1
         vn = vn + rvn
         if (abs(rvn/(vn + 1.0e-50_dp)) < eps4 .and. j > ia/2) exit
      end do
      vn = 2.0_dp*vn
      if (iff == 1) vn = vn*1.0e300_dp
   end subroutine cdcof4a

!=============================================================================
!  cdrb12cc -- R^(1), R^(2) as series of spherical Bessel functions
!=============================================================================
   subroutine cdrb12cc(r, dr, ar, adr, id, m, n, c, x, ad, ia, ier)
      complex(dp), intent(out) :: r, dr, ar, adr
      integer, intent(in) :: id, m, n, ia
      complex(dp), intent(in) :: c, ad(ia)
      real(dp), intent(in) :: x
      integer, intent(out) :: ier

      real(dp), parameter :: factor = 1.0e250_dp
      complex(dp) :: besj(4*nterms+1), besy(4*nterms+1)
      integer :: iiy(4*nterms+1)
      complex(dp) :: ac, sig, rk, f2, rs, rd, br, abr
      real(dp) :: am, pi2, x2, f1
      integer :: m2, jk, i, j, i1, num, nm2
      logical :: skip_ar

      r = czero; dr = czero; ar = czero; adr = czero
      ier = 0
      if (x < 1.0e-10_dp) then
         ier = 10
         return
      end if
      m2 = 2*m
      if (mod(n - m, 2) /= 0) then
         jk = 1
      else
         jk = 0
      end if
      am = m
      am = am/2.0_dp
      pi2 = pi/2.0_dp
      x2 = x*x
      ac = c*x
      sig = czero
      do i = 1, ia
         j = 2*(i - 1) + jk
         f2 = ad(i)
         do i1 = 1, m2
            f2 = f2*(j + i1)
         end do
         sig = sig + f2
      end do
      rk = sqrt(pi2/ac)/sig*((x2 - ak)/x2)**am
      num = int(3*abs(ac) + 10)
      if (x > 2.0_dp) num = int(3*abs(ac) + 6*(m + x))
      if (num < 60) num = 60
      if (num - 1 >= 4*nterms) write (*, *) num, 4*nterms
      if (num - 1 >= 4*nterms) num = 4*nterms

      if (id == 1) call cessel0(1.0_dp/ac, num, besj, besy, iiy)
      if (id == 0) call bessjj0(1.0_dp/ac, num, besj)

      f1 = ak/x*m/(x2 - ak)
      nm2 = min(ia, num/2 - 2)
      skip_ar = .false.
      rs = czero; rd = czero; br = czero; abr = czero

      do i = 1, nm2
         j = 2*(i - 1) + jk
         f2 = ad(i)*ci**(j + m - n)
         do i1 = 1, m2
            f2 = f2*(j + i1)
         end do
         if (i /= 1) then
            if (abs(rd) < eps3 .and. abs(rs) < eps3 .and. i > (m + 5) .and. id == 0) then
               skip_ar = .true.
               exit
            end if
         end if
         rs = f2*besj(j+m+1)
         rd = f2*((f1 + (m + j)/x)*besj(j+m+1) - c*besj(j+m+2))
         r = r + rs
         dr = dr + rd
         if (id == 0) cycle
         if (i /= 1) then
            if (abs(br) < eps3 .and. abs(abr) < eps3 .and. i > (m + 5)) exit
         end if
         br = besy(j+m+1)*f2*factor**(iiy(j+m+1))
         ar = ar + br
         abr = (f1 + (m + j)/x)*besy(j+m+1)*f2*factor**(iiy(j+m+1)) &
               - c*besy(j+m+2)*f2*factor**(iiy(j+m+1))
         adr = adr + abr
      end do
      if (.not. skip_ar) then
         ar = ar*rk
         adr = adr*rk
      end if
      r = r*rk
      dr = dr*rk
   end subroutine cdrb12cc

!=============================================================================
!  cdrf12cc -- R^(1), R^(2) as series of associated Legendre functions
!=============================================================================
   subroutine cdrf12cc(r, dr, ar, adr, i12, m, n, c, ad, ia, bd, ib, ier)
      complex(dp), intent(out) :: r, dr
      complex(dp), intent(inout) :: ar, adr
      integer, intent(in) :: i12, m, n, ia, ib
      complex(dp), intent(in) :: c, ad(ia), bd(*)
      integer, intent(out) :: ier

      complex(dp) :: sig, rf2, rk1, rk2, rs, rd, rbr, rabr
      real(dp) :: eps4, am2, ff3
      integer :: m2, nm1, nm2, jk, jn, jm, i, j, j1, j2, mr

      eps4 = 1.0e-20_dp
      ier = 0
      m2 = 2*m
      am2 = 2.0_dp*m
      nm1 = n - m
      nm2 = n + m
      jk = merge(1, 0, mod(n - m, 2) /= 0)
      jn = merge(1, 0, mod(n, 2) /= 0)
      jm = merge(1, 0, mod(m, 2) /= 0)

      sig = czero
      do i = 1, ia
         j = 2*(i - 1) + jk
         rf2 = ad(i)*fact2(m2 + j, j)
         sig = sig + rf2
         if (abs(rf2/sig) < eps4) exit
      end do
      if ((nm1 - jk)/2 + 1 > 170) error stop 1702
      ff3 = fact3(nm2 + jk, (nm2 + jk)/2, 1.0_dp/fact((nm1 - jk)/2 + 1))
      rk1 = (am2 + 1.0_dp + jk*2.0_dp)/2.0_dp**nm2*ff3/ &
            ad(1)/fact(m+1)*sig/c**(m + jk)
      if (kob == 1) rk1 = rk1/(-ci)**(m + jk)*(-ci)**jn
      r = czero
      dr = czero

      if (i12 /= 0) then
         rk2 = 2.0_dp**nm1/(am2 - 1.0_dp)*fact(m2+1)/fact(m+1)/ff3*bd(m+1)*sig/c**(m - 1 - jk)
         if (jk == 1 .and. kob == 0) rk2 = -rk2/(am2 - 3.0_dp)
         if (kob /= 0) then
            if (jk == 1 .and. jm == 0) rk2 = rk2/(am2 - 3.0_dp)
            if (jk == 1 .and. jm == 1) rk2 = -rk2/(am2 - 3.0_dp)*(-ci)
            if (jk == 0 .and. jm == 0) rk2 = rk2*ci
            if (kob == 1) rk2 = rk2/(-ci)**(m - 1 - jk)
         end if
         ar = czero
         adr = czero
         do i = 1, m
            j2 = 2*(m - i) + 1 + jk
            adr = adr + bd(i+1)*qd_leg(j2)
            ar = ar + bd(i+1)*q_leg(j2)
         end do
      end if

      do i = 1, ia
         j1 = 2*(i - 1) + 1 + jk
         j2 = 2*(m + i) - 1 + jk
         rs = ad(i)*p_leg(j1)
         rd = ad(i)*pd_leg(j1)
         r = r + rs
         dr = dr + rd
         if (i12 /= 0) then
            rbr = ad(i)*q_leg(j2)
            rabr = ad(i)*qd_leg(j2)
            ar = ar + rbr
            adr = adr + rabr
            if (abs(rbr/ar) < eps4 .and. abs(rabr/adr) < eps3) exit
         end if
         if (abs(rs/r) < eps4 .and. abs(rd/dr) < eps4) exit
      end do

      if (i12 /= 0) then
         mr = m + 1
         do i = mr, ib
            j2 = 2*(i - m) - jk
            if (kob == 1 .and. mod(m, 2) == 0) then
               rbr = -bd(i+1)*p_leg(j2)
               rabr = -bd(i+1)*pd_leg(j2)
            else
               rbr = bd(i+1)*p_leg(j2)
               rabr = bd(i+1)*pd_leg(j2)
            end if
            ar = ar + rbr
            adr = adr + rabr
            if (abs(rbr/ar) < eps4 .and. abs(rabr/adr) < eps4) exit
         end do
         ar = ar/rk2
         adr = adr/rk2
      end if
      r = r/rk1
      dr = dr/rk1
   end subroutine cdrf12cc

!=============================================================================
!  cdrg1cn -- prolate R^(1), R^(2) via the Jaffe-type expansion
!=============================================================================
   subroutine cdrg1cn(r, dr, ar, adr, i12, m, n, c, x, ad, ia, bd, ib, rrl, ier)
      complex(dp), intent(out) :: r, dr
      complex(dp), intent(inout) :: ar, adr
      integer, intent(in) :: i12, m, n, ia, ib
      complex(dp), intent(in) :: c, ad(ia), bd(*), rrl
      real(dp), intent(in) :: x
      integer, intent(out) :: ier

      complex(dp) :: sg(1000), sa(1000), sb(1000), sam(100)
      complex(dp) :: rl, rp, re, rp2, rrp2l, rg0, sig, rf2, rk2
      complex(dp) :: sssa0, sssa1, rf, rf1, rf3, rf4, rf11, rs2, rrs2, rd21, rd23
      complex(dp) :: sg0, drr, drr2, drr3, r2, rr2, dr2, s2ic, sa00, sb0, rb, rs, rd
      real(dp) :: aksi, eps4, am, xx0, xx, ff3, ff2, a, g
      integer :: m2, nm1, nm2, jk, jm, i, j, jj, jjj, ir, irm, mr, ii

      ier = 0
      rk2 = czero
      ! lambda according to Slavyanov
      rl = rrl - c**2
      aksi = 0.0_dp
      if (kob == 0) aksi = x**2 - 1.0_dp
      if (kob == 1) aksi = x**2 + 1.0_dp
      eps4 = 1.0e-20_dp
      am = 0.5_dp*m
      xx0 = 0.0_dp
      xx = 0.0_dp
      if (kob == 0) xx0 = ((x - 1.0_dp)/(x + 1.0_dp))**am
      if (kob == 0) xx = (x - 1.0_dp)/(x + 1.0_dp)

      rp = -ci*c
      re = exp(-rp*(x - 1.0_dp))
      rp2 = 2.0_dp*rp
      rrp2l = (rp2 + 1.0_dp)*(m + 1.0_dp) + rl

      ! coefficients
      m2 = 2*m
      nm1 = n - m
      nm2 = n + m
      jk = merge(1, 0, mod(n - m, 2) /= 0)

      ff3 = fact3(nm2 + jk, (nm2 + jk)/2, 1.0_dp/fact((nm1 - jk)/2 + 1))
      rg0 = (2.0_dp*m + 1.0_dp + jk*2.0_dp)/2.0_dp**(nm2 + 1.0_dp)*ff3/ad(1)/c**(m + jk)
      rg0 = 1.0_dp/rg0

      if (i12 /= 0) then
         jm = merge(1, 0, mod(m, 2) /= 0)
         sig = czero
         do i = 1, ia
            j = 2*(i - 1) + jk
            rf2 = ad(i)*fact2(m2 + j, j)
            sig = sig + rf2
            if (abs(rf2/sig) < eps4) exit
         end do
         rk2 = 2.0_dp**nm1/(m2 - 1.0_dp)*fact(m2+1)/fact(m+1)/ff3*bd(m+1)*sig/c**(m - 1 - jk)
         if (jk == 1 .and. kob == 0) rk2 = -rk2/(m2 - 3.0_dp)
         if (kob /= 0) then
            if (jk == 1 .and. jm == 0) rk2 = rk2/(m2 - 3.0_dp)
            if (jk == 1 .and. jm == 1) rk2 = -rk2/(m2 - 3.0_dp)*(-ci)
            if (jk == 0 .and. jm == 0) rk2 = rk2*ci
            if (kob == 1) rk2 = rk2/(-ci)**(m - 1 - jk)
         end if
      end if

      ! ---- 1s ----
      sssa0 = czero
      do i = 1, m
         j = -2*(m - i + 1) + jk
         if (j >= -m) exit
         sssa0 = sssa0 + bd(m-i+2)*(-1)**(-j-1)*fact(-j)*fact(m2+j+1)
      end do

      ! ---- 2s ----
      rf1 = czero
      do i = 1, ia
         ir = 2*(i - 1) + jk
         rf3 = czero
         do j = 1, m
            ff2 = fact(m+1)/fact(m-j+1)
            rf2 = czero
            do jj = 1, j + 1
               jjj = jj - 1
               rf2 = rf2 + fact2(m, m - jjj)**2/fact(jjj+1)/fact(j-jjj+1)*fact2(ir + m2 - jjj, ir + jjj)
            end do
            rf2 = rf2*(-1)**j/(2.0_dp*j)*ff2
            rf3 = rf3 + (-1)**(j - 1)/(2.0_dp*j)*ff2**2*fact2(ir + m2 - j, ir + j) + rf2
         end do
         irm = (ir + m - 1 + jk)/2
         rf4 = czero
         if (ir /= 0) then
            do j = 1, irm
               jj = j - 1
               if (2*jj > ir - 1) exit
               rf4 = rf4 + (2.0_dp*ir + m2 - 4.0_dp*jj - 1.0_dp)/(ir + m - jj) &
                     /(2.0_dp*jj + 1.0_dp)*fact2(ir + m2 - 2*jj - 1, ir - 2*jj - 1)
            end do
         end if
         rf11 = ad(i)*(rf3 - rf4)
         rf1 = rf1 + rf11
         if (abs(rf11/rf1) < eps4) exit
      end do
      sssa0 = sssa0 + rf1

      ! ---- 3s ----
      rf1 = czero
      mr = m + 1
      do i = mr, ib
         j = 2*i - jk
         rf11 = bd(i+1)*fact2(j - 1, j - m2 - 1)
         rf1 = rf1 + rf11
         if (abs(rf11/rf1) < eps4) exit
      end do
      sssa0 = sssa0 + rf1

      ! ---- function of the first kind: start of the recurrence ----
      sg0 = (1.0_dp, 0.0_dp)
      sg(1) = rrp2l/(m + 1.0_dp)*sg0
      sg(2) = ((2.0_dp*(m + 2.0_dp + rp2) + rrp2l)*sg(1) - (m + 1.0_dp)*sg0)/(2.0_dp*(m + 2.0_dp))
      r = sg0 + sg(1)*xx + sg(2)*xx**2
      drr = sg(1)*1.0_dp + sg(2)*2.0_dp*xx

      ! ---- function of the second kind: start of the recurrence ----
      r2 = czero; rr2 = czero; drr2 = czero; drr3 = czero
      if (i12 /= 0) then
         sam = czero
         sam(1) = -2.0_dp/m/c/rg0
         if (-m + 1 /= 0) then
            sam(2) = (2.0_dp*m*(rp2 + 1.0_dp) - rrp2l)/(m - 1.0_dp)*sam(1)
            if (-m + 2 /= 0) then
               do i = 1, m - 2
                  j = -m + i
                  sam(i+2) = ((2.0_dp*j*(j + rp2 + m + 1.0_dp) + rrp2l)*sam(i+1) - &
                              (j*(j + m)*sam(i)))/((j + 1.0_dp)*(j + m + 1.0_dp))
               end do
            end if
         end if

         ! a0
         s2ic = 2.0_dp*ci*c
         sssa1 = czero
         do i = 1, m
            j = i - 1
            rf = czero
            do ii = 1, m - j
               rf = rf + s2ic**ii/fact(ii+1)*fact(m-j)/fact(ii)/fact(m-ii-j+1)
            end do
            sssa1 = sssa1 + sam(i)*rf
         end do

         sa00 = 2.0_dp/fact(m+1)*sssa0/rk2 - sssa1

         sb0 = (-2.0_dp*(rp2 + m) + rrp2l)*sam(m)/m
         if (m > 1) sb0 = sb0 + (m - 1.0_dp)*sam(m-1)/m

         sb(1) = rrp2l/(m + 1.0_dp)*sb0
         sb(2) = ((2.0_dp*(m + 2.0_dp + rp2) + rrp2l)*sb(1) - (m + 1.0_dp)*sb0)/(2.0_dp*(m + 2.0_dp))
         sa(1) = (rrp2l*sa00 - (m + 2.0_dp)*sb(1) + 2.0_dp*(rp2 + m + 1.0_dp)*sb0)/(m + 1.0_dp)
         sa(2) = ((2.0_dp*(m + 2.0_dp + rp2) + rrp2l)*sa(1) - (m + 1.0_dp)*sa00 - ((m + 4.0_dp)*sb(2) &
                 - 2.0_dp*(rp2 + m + 3.0_dp)*sb(1) + (m + 2.0_dp)*sb0))/(2.0_dp*(m + 2.0_dp))

         do i = 1, m
            r2 = r2 + sam(m-i+1)/xx**i
            drr2 = drr2 - i*sam(m-i+1)/xx**(i + 1)
         end do

         r2 = r2 + sa00 + sa(1)*xx + sa(2)*xx**2
         rr2 = sb0 + sb(1)*xx + sb(2)*xx**2
         drr2 = drr2 + sa(1)*1.0_dp + 2.0_dp*sa(2)*xx
         drr3 = sb(1)*1.0_dp + 2.0_dp*sb(2)*xx
      end if

      ! ---- main recurrence ----
      do i = 2, 998
         a = (i + 1.0_dp)*(i + m + 1.0_dp)
         rb = 2.0_dp*i*(i + m + 1.0_dp + rp2) + rrp2l
         g = i*(i + m)

         sg(i+1) = (rb*sg(i) - g*sg(i-1))/a
         rs = sg(i+1)*xx**(i + 1.0_dp)
         r = r + rs
         rd = sg(i+1)*(i + 1.0_dp)*xx**i
         drr = drr + rd

         rs2 = czero; rrs2 = czero; rd21 = czero; rd23 = czero
         if (i12 /= 0) then
            if (abs(sb(i)) > 1.0e300_dp) return
            sb(i+1) = (rb*sb(i) - g*sb(i-1))/a
            sa(i+1) = (rb*sa(i) - g*sa(i-1) - ((2.0_dp*i + m + 2.0_dp)* &
                      sb(i+1) - 2.0_dp*(2.0_dp*i + rp2 + m + 1.0_dp)*sb(i) + &
                      (2.0_dp*i + m)*sb(i-1)))/a
            rs2 = sa(i+1)*xx**(i + 1.0_dp)
            rrs2 = sb(i+1)*xx**(i + 1.0_dp)
            r2 = r2 + rs2
            rr2 = rr2 + rrs2
            rd21 = sa(i+1)*(i + 1.0_dp)*xx**i
            rd23 = sb(i+1)*(i + 1.0_dp)*xx**i
            drr2 = drr2 + rd21
            drr3 = drr3 + rd23
         end if
         if (abs(rs/r) < eps4 .and. abs(rd/drr) < eps4) exit
      end do

      r = r*xx0/(x + 1.0_dp)*re*rg0
      dr = r*(-rp - 1.0_dp/(x + 1.0_dp) + m/aksi) + &
           drr*(2.0_dp/(x + 1.0_dp)**3)*xx0*re*rg0

      if (i12 /= 0) then
         r2 = (rr2*log(xx) + r2)*xx0/(x + 1.0_dp)*re
         dr2 = r2*(-rp - 1.0_dp/(x + 1.0_dp) + m/aksi) + &
               (drr2 + rr2/xx + drr3*log(xx))*(2.0_dp/(x + 1.0_dp)**3)*xx0*re
         ar = r2
         adr = dr2
      end if

      if (aimag(c) == 0.0_dp) then
         r = cmplx(real(r, dp), 0.0_dp, dp)
         dr = cmplx(real(dr, dp), 0.0_dp, dp)
         ar = cmplx(real(ar, dp), 0.0_dp, dp)
         adr = cmplx(real(adr, dp), 0.0_dp, dp)
      end if
   end subroutine cdrg1cn

!=============================================================================
!  Runge-Kutta-Fehlberg (4,5) integrator for complex systems (NEQN <= 4)
!=============================================================================
   subroutine cdrkf45a(f, neqn, y, t, tout, relerr, abserr, iflag)
      procedure(rhs_fn) :: f
      integer, intent(in) :: neqn
      complex(dp), intent(inout) :: y(neqn)
      real(dp), intent(inout) :: t, relerr, abserr
      real(dp), intent(in) :: tout
      integer, intent(inout) :: iflag
      call cdrkfsa(f, neqn, y, t, tout, relerr, abserr, iflag)
   end subroutine cdrkf45a

   ! The integrator keeps its state between calls (as all RKF45 codes do),
   ! hence the SAVE statement.
   subroutine cdrkfsa(f, neqn, y, t, tout, relerr, abserr, iflag)
      procedure(rhs_fn) :: f
      integer, intent(in) :: neqn
      complex(dp), intent(inout) :: y(neqn)
      real(dp), intent(inout) :: t, relerr, abserr
      real(dp), intent(in) :: tout
      integer, intent(inout) :: iflag

      real(dp), parameter :: remin = 1.0e-16_dp
      integer, parameter :: maxnfe = 500000
      logical :: hfaild, output
      complex(dp) :: yp(4), f1(4), f2(4), f3(4), f4(4), f5(4), fout(4)
      real(dp) :: h, savre, savae, a, ae, dt, ee, eeoet, esttol, et, hmin
      real(dp) :: rer, s, scl, tol, toln, u26, eps, ypk
      integer :: mflag, kflag, init, kop, nfe, jflag, k
      save

      if (neqn < 1) goto 10
      if (relerr < 0.0_dp .or. abserr < 0.0_dp) goto 10
      mflag = abs(iflag)
      if (mflag == 0 .or. mflag > 8) goto 10
      if (mflag /= 1) goto 20
      ! unit roundoff (the original found it with a halving loop)
      eps = epsilon(1.0_dp)/2.0_dp
      u26 = 26.0_dp*eps
      goto 50
10    iflag = 8
      return
20    if (t == tout .and. kflag /= 3) goto 10
      if (mflag /= 2) goto 25
      if (kflag == 3 .or. init == 0) goto 45
      if (kflag == 4) goto 40
      if (kflag == 5 .and. abserr == 0.0_dp) goto 30
      if (kflag == 6 .and. relerr <= savre .and. abserr <= savae) goto 30
      goto 50
25    if (iflag == 3) goto 45
      if (iflag == 4) goto 40
      if (iflag == 5 .and. abserr > 0.0_dp) goto 45
30    error stop 'cdrkfsa: integration cannot be continued'
40    nfe = 0
      if (mflag == 2) goto 50
45    iflag = jflag
      if (kflag == 3) mflag = abs(iflag)
50    jflag = iflag
      kflag = 0
      savre = relerr
      savae = abserr
      rer = 2.0_dp*eps + remin
      if (relerr >= rer) goto 55
      relerr = rer
      iflag = 3
      kflag = 3
      return
55    dt = tout - t
      if (mflag == 1) goto 60
      if (init == 0) goto 65
      goto 80
60    init = 0
      kop = 0
      a = t
      call f(a, y, yp)
      nfe = 1
      if (t /= tout) goto 65
      iflag = 2
      return
65    init = 1
      h = abs(dt)
      toln = 0.0_dp
      do k = 1, neqn
         tol = relerr*abs(y(k)) + abserr
         if (tol <= 0.0_dp) cycle
         toln = tol
         ypk = abs(yp(k))
         if (ypk*h**5 > tol) h = (tol/ypk)**0.2_dp
      end do
      if (toln <= 0.0_dp) h = 0.0_dp
      h = 1.0e-4_dp
      h = max(h, u26*max(abs(t), abs(dt)))
      jflag = sign(2, iflag)
80    h = sign(h, dt)
      if (abs(h) >= 2.0_dp*abs(dt)) kop = kop + 1
      if (kop /= 100) goto 85
      kop = 0
      iflag = 7
      return
85    if (abs(dt) > u26*abs(t)) goto 95
      do k = 1, neqn
         y(k) = y(k) + dt*yp(k)
      end do
      a = tout
      call f(a, y, yp)
      nfe = nfe + 1
      goto 300
95    output = .false.
      scl = 2.0_dp/relerr
      ae = scl*abserr
100   hfaild = .false.
      hmin = u26*abs(t)
      dt = tout - t
      if (abs(dt) >= 2.0_dp*abs(h)) goto 200
      if (abs(dt) > abs(h)) goto 150
      output = .true.
      h = dt
      goto 200
150   h = 0.5_dp*dt
200   if (nfe <= maxnfe) goto 220
      iflag = 4
      kflag = 4
      return
220   call cfehl(f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, fout)
      f1(1:neqn) = fout(1:neqn)
      nfe = nfe + 5
      eeoet = 0.0_dp
      do k = 1, neqn
         et = abs(y(k)) + abs(f1(k)) + ae
         if (et <= 0.0_dp) then
            iflag = 5
            return
         end if
         ee = abs((-2.09e3_dp*yp(k) + (2.197e4_dp*f3(k) - &
              1.5048e4_dp*f4(k))) + (2.2528e4_dp*f2(k) - 2.736e4_dp*f5(k)))
         eeoet = max(eeoet, ee/et)
      end do
      esttol = abs(h)*eeoet*scl/7.524e5_dp
      if (esttol <= 1.0_dp) goto 260
      hfaild = .true.
      output = .false.
      s = 0.1_dp
      if (esttol < 5.9049e4_dp) s = 0.9_dp/esttol**0.2_dp
      h = s*h
      if (abs(h) > hmin) goto 200
      iflag = 6
      kflag = 6
      return
260   t = t + h
      do k = 1, neqn
         y(k) = f1(k)
      end do
      a = t
      call f(a, y, yp)
      nfe = nfe + 1
      s = 5.0_dp
      if (esttol > 1.889568e-4_dp) s = 0.9_dp/esttol**0.2_dp
      if (hfaild) s = min(s, 1.0_dp)
      h = sign(max(s*abs(h), hmin), h)
      if (output) goto 300
      if (iflag > 0) goto 100
      iflag = -2
      return
300   t = tout
      iflag = 2
   end subroutine cdrkfsa

   ! One Fehlberg step; result returned in s(:)
   subroutine cfehl(f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s)
      procedure(rhs_fn) :: f
      integer, intent(in) :: neqn
      complex(dp), intent(in) :: y(neqn), yp(neqn)
      real(dp), intent(in) :: t, h
      complex(dp), intent(inout) :: f1(neqn), f2(neqn), f3(neqn), f4(neqn), f5(neqn)
      complex(dp), intent(out) :: s(neqn)
      real(dp) :: ch
      integer :: k

      ch = h/4.0_dp
      do k = 1, neqn
         f5(k) = y(k) + ch*yp(k)
      end do
      call f(t + ch, f5, f1)
      ch = 3.0_dp*h/3.2e1_dp
      do k = 1, neqn
         f5(k) = y(k) + ch*(yp(k) + 3.0_dp*f1(k))
      end do
      call f(t + 3.0_dp*h/8.0_dp, f5, f2)
      ch = h/2.197e3_dp
      do k = 1, neqn
         f5(k) = y(k) + ch*(1.932e3_dp*yp(k) + (7.296e3_dp*f2(k) - 7.2e3_dp*f1(k)))
      end do
      call f(t + 1.2e1_dp*h/1.3e1_dp, f5, f3)
      ch = h/4.104e3_dp
      do k = 1, neqn
         f5(k) = y(k) + ch*((8.341e3_dp*yp(k) - 8.45e2_dp*f3(k)) &
                 + (2.944e4_dp*f2(k) - 3.2832e4_dp*f1(k)))
      end do
      call f(t + h, f5, f4)
      ch = h/2.052e4_dp
      do k = 1, neqn
         f1(k) = y(k) + ch*((-6.08e3_dp*yp(k) + (9.295e3_dp*f3(k) - &
                 5.643e3_dp*f4(k))) + (4.104e4_dp*f1(k) - 2.8352e4_dp*f2(k)))
      end do
      call f(t + h/2.0_dp, f1, f5)
      ch = h/7.61805e6_dp
      do k = 1, neqn
         s(k) = y(k) + ch*((9.0288e5_dp*yp(k) + (3.855735e6_dp*f3(k) - &
                1.371249e6_dp*f4(k))) + (3.953664e6_dp*f2(k) + 2.7702e5_dp*f5(k)))
      end do
   end subroutine cfehl

!=============================================================================
!  cdrsf20 -- prolate R^(2) by integrating the radial equation (Prufer-type
!             variables) from large xi down to x
!=============================================================================
   subroutine cdrsf20(r, dr, x, r1, r2, m, c, rl, ier)
      complex(dp), intent(out) :: r, dr
      real(dp), intent(in) :: x
      complex(dp), intent(in) :: r1, r2, c, rl
      integer, intent(in) :: m
      integer, intent(out) :: ier

      character(*), parameter :: fmt31 = '(" TOLERANCES RESET",2d12.3)'
      complex(dp) :: y(4), rra, rrb, r3
      real(dp) :: pl, ac, tksi, t, elerr, abserr, tfinal, tout
      integer :: neqn, iflag

      f2_m = m
      f2_c = c
      f2_rl = rl
      neqn = 4
      pl = real(rl, dp)
      ac = real(c, dp)

      tksi = sqrt((1.0_dp + pl/ac**2)/2.0_dp + sqrt((1.0_dp - pl/ac**2)**2/4.0_dp + &
             (m**2 - 1.0_dp)/ac**2))
      if (m == 1) tksi = 1.0_dp
      t = tksi + x/ac*10.0_dp
      elerr = 1.0e-8_dp
      abserr = 0.0_dp
      tfinal = x
      iflag = 1
      tout = t
      ier = 0
      call cdrfp0(m, c, rl, t, y)

      do
         call cdrkf45a(cdrfp2, neqn, y, t, tout, elerr, abserr, iflag)
         rra = y(1) + y(2)
         rrb = y(3) + y(4)
         r3 = -ci/(c*aksi_g*((rrb/rra + x/aksi_g)*r1 + r2))
         r = -ci*r3
         dr = -(x/aksi_g + rrb/rra)*r + r1*ci*rrb/rra

         select case (iflag)
         case (2)
            tout = x
            if (t > tfinal) cycle
            ier = 0
            exit
         case (3)
            write (lu, fmt31) elerr, abserr
            ier = 3
         case (4)
            write (lu, '(a)') ' MANY STEPS'
            ier = 4
         case (5)
            abserr = 1.0e-9_dp
            write (lu, fmt31) elerr, abserr
            ier = 5
         case (6)
            elerr = 10.0_dp*elerr
            write (lu, fmt31) elerr, abserr
            iflag = 2
            ier = 6
         case (7)
            write (lu, '(a)') ' MUCH OUTPUT'
            iflag = 2
            ier = 7
         case default
            write (lu, '(a)') ' IMPROPER CALL'
            ier = 8
            exit
         end select
      end do
   end subroutine cdrsf20

   ! initial values for cdrsf20
   subroutine cdrfp0(m, c, rl, ksi, y)
      integer, intent(in) :: m
      complex(dp), intent(in) :: c, rl
      real(dp), intent(in) :: ksi
      complex(dp), intent(out) :: y(4)
      complex(dp) :: c2, gksi, al, an
      real(dp) :: am2, aksi

      am2 = m*m
      c2 = c**2
      aksi = ksi**2 - 1.0_dp
      gksi = c2 - (rl - c2)/aksi + (1.0_dp - am2)/aksi**2
      al = ksi*(2.0_dp*(1.0_dp - am2) - (rl - c2)*aksi)/(2.0_dp*gksi*aksi**3)
      an = sqrt(1.0_dp + abs(gksi) + al**2)
      y(1) = 1.0_dp/an
      y(2) = czero
      y(3) = -al/an
      y(4) = -sqrt(abs(gksi))/an
   end subroutine cdrfp0

   ! right-hand side for cdrsf20
   subroutine cdrfp2(x, y, yp)
      real(dp), intent(in) :: x
      complex(dp), intent(in) :: y(*)
      complex(dp), intent(out) :: yp(*)
      complex(dp) :: c2, gksi, akk
      real(dp) :: aksi

      aksi = x**2 - 1.0_dp
      c2 = f2_c**2
      gksi = c2 - (f2_rl - c2)/aksi + (1.0_dp - f2_m**2)/aksi**2
      ! (operator precedence exactly as in the original)
      akk = (y(1)*y(3) + y(2)*y(4)*(1.0_dp - gksi)/(y(1)**2 + y(2)**2 + y(3)**2 + &
             y(4)**2))
      yp(1) = akk*y(1) - y(3)
      yp(2) = akk*y(2) - y(4)
      yp(3) = akk*y(3) + gksi*y(1)
      yp(4) = akk*y(4) + gksi*y(2)
   end subroutine cdrfp2

!=============================================================================
!  cdrsf212 -- oblate R^(1), R^(2) by integrating the radial equation from 0
!=============================================================================
   subroutine cdrsf212(r, dr, ar, adr, i12, m, n, c, x, rl, ad, ia, bd, ib, ier)
      complex(dp), intent(inout) :: r, dr, ar, adr
      integer, intent(in) :: i12, m, n, ia, ib
      complex(dp), intent(in) :: c, rl, ad(ia), bd(ib)
      real(dp), intent(in) :: x
      integer, intent(out) :: ier

      character(*), parameter :: fmt31 = '(" TOLERANCES RESET",2d12.3)'
      complex(dp) :: y(4)
      real(dp) :: t, elerr, abserr, tfinal, tout
      integer :: neqn, iflag

      f2_m = m
      f2_c = c
      f2_rl = rl
      f2_i12 = i12

      neqn = 2
      if (i12 == 2) neqn = 4
      t = 0.0_dp
      elerr = 1.0e-15_dp
      abserr = 0.0_dp
      tfinal = x
      iflag = 1
      tout = t
      ier = 0
      call cdrf2k02(m, n, c, ad, ia, bd, ib, i12, y)

      do
         call cdrkf45a(cdrf222, neqn, y, t, tout, elerr, abserr, iflag)

         if (i12 == 0) then
            r = y(1)
            dr = y(2)
         end if
         if (i12 == 1) then
            ar = y(1)
            adr = y(2)
         end if
         if (i12 == 2) then
            r = y(1)
            dr = y(2)
            ar = y(3)
            adr = y(4)
         end if

         select case (iflag)
         case (2)
            tout = x
            if (t < tfinal) cycle
            ier = 0
            exit
         case (3)
            write (lu, fmt31) elerr, abserr
            ier = 3
         case (4)
            write (lu, '(a)') ' MANY STEPS'
            ier = 4
         case (5)
            abserr = 1.0e-9_dp
            write (lu, fmt31) elerr, abserr
            ier = 5
         case (6)
            elerr = 10.0_dp*elerr
            write (lu, fmt31) elerr, abserr
            iflag = 2
            ier = 6
         case (7)
            write (lu, '(a)') ' MUCH OUTPUT'
            iflag = 2
            ier = 7
         case default
            write (lu, '(a)') ' IMPROPER CALL'
            ier = 8
            exit
         end select
      end do
   end subroutine cdrsf212

   ! right-hand side for cdrsf212
   subroutine cdrf222(x, y, yp)
      real(dp), intent(in) :: x
      complex(dp), intent(in) :: y(*)
      complex(dp), intent(out) :: yp(*)
      complex(dp) :: ar
      real(dp) :: aksi, am

      aksi = x**2 + 1.0_dp
      am = f2_m*f2_m
      ar = f2_rl - (f2_c*x)**2 - am/aksi

      if (f2_i12 == 0 .or. f2_i12 == 1) then
         yp(1) = y(2)
         yp(2) = (-2.0_dp*x*y(2) + ar*y(1))/aksi
      end if
      if (f2_i12 == 2) then
         yp(1) = y(2)
         yp(2) = (-2.0_dp*x*y(2) + ar*y(1))/aksi
         yp(3) = y(4)
         yp(4) = (-2.0_dp*x*y(4) + ar*y(3))/aksi
      end if
   end subroutine cdrf222

   ! initial values (at xi = 0) for cdrsf212
   subroutine cdrf2k02(m, n, c, ad, ia, bd, ib, i12, y)
      integer, intent(in) :: m, n, ia, ib, i12
      complex(dp), intent(in) :: c, ad(ia), bd(ib)
      complex(dp), intent(out) :: y(4)

      complex(dp) :: r, r1, r2, s2, t, ss, ss2
      real(dp) :: f, eps4
      integer :: m2, jk, i, j, jj, mf, jmn, jmn2

      eps4 = 1.0e-20_dp
      m2 = 2*m
      jk = merge(1, 0, mod(n - m, 2) /= 0)
      s2 = czero
      do i = 1, ia
         j = 2*(i - 1) + jk
         f = 1.0_dp
         do jj = 1, m2
            f = f*(j + jj)
         end do
         t = ad(i)*f
         s2 = s2 + t
         if (abs(t/s2) < eps4) exit
      end do
      mf = n + m + 1 + jk
      f = fact(mf)

      ! phase factors (computed always; the original left ss undefined for i12 = 0)
      jmn = mod(n - m - jk, 4)
      ss = (1.0_dp, 0.0_dp)
      if (jmn == 1) ss = ci
      if (jmn == 2) ss = (-1.0_dp, 0.0_dp)
      if (jmn == 3) ss = -ci
      jmn2 = abs(mod(m - n - jk, 4))
      ss2 = (1.0_dp, 0.0_dp)
      if (jmn2 == 1) ss2 = ci
      if (jmn2 == 2) ss2 = (-1.0_dp, 0.0_dp)
      if (jmn2 == 3) ss2 = -ci

      r = czero; r1 = czero; r2 = czero
      if (i12 == 1 .or. i12 == 2) then
         r = ss*(m2 - 1.0_dp)*fact(m+1)*pi/2.0_dp**(2*n - m + 1)*c**(m - 1 - jk)/ &
             fact(m2+1)/bd(m+1)/s2* &
             (f/fact((n - m - jk)/2 + 1)/fact((n + m + jk)/2 + 1))**2
         r1 = ss2*s2/2.0_dp**m/c**(m + 1 + jk)/ad(1)/fact(m+1)
      end if
      if (i12 == 0 .or. i12 == 2) then
         r2 = ss*2.0_dp**m*c**(m + jk)*ad(1)*fact(m+1)/s2
      end if

      y = czero
      if (jk /= 1) then
         if (i12 == 0) then
            y(1) = r2/(m2 + 1.0_dp)
            y(2) = czero
         end if
         if (i12 == 1) then
            y(1) = r
            y(2) = r1*(m2 + 1.0_dp)
         end if
         if (i12 == 2) then
            y(1) = r2/(m2 + 1.0_dp)
            y(2) = czero
            y(3) = r
            y(4) = r1*(m2 + 1.0_dp)
         end if
      else
         if (i12 == 0) then
            y(1) = czero
            y(2) = r2/(m2 + 3.0_dp)
         end if
         if (i12 == 1) then
            y(1) = r1*(m2 + 3.0_dp)
            y(2) = r*(m2 - 3.0_dp)
         end if
         if (i12 == 2) then
            y(1) = czero
            y(2) = r2/(m2 + 3.0_dp)
            y(3) = r1*(m2 + 3.0_dp)
            y(4) = r*(m2 - 3.0_dp)
         end if
      end if
   end subroutine cdrf2k02

!=============================================================================
!  Helpers shared by the four eigenvalue routines
!=============================================================================

   ! Continued-fraction residual SL(lambda) = RAA - RBB (the 'from below to
   ! top' and 'from top to below' iterations).  guard = .true. adds 1e-50 to
   ! the denominator, as was done in CDLAMo_l only.
   function cf_residual(lam, m, n1, nm2, nmax, rc2, guard) result(sl)
      complex(dp), intent(in) :: lam, rc2
      integer, intent(in) :: m, n1, nm2, nmax
      logical, intent(in) :: guard
      complex(dp) :: sl

      complex(dp) :: rda(max(nmax, 3)), raa, den
      real(dp) :: am, amm, ajj, aii, aj
      integer :: i, i1, j, m2

      m2 = m*m
      am = m
      amm = 2.0_dp*m

      ! iteration 'from below to top'
      if (mod(n1, 2) == 0) then
         rda(2) = (amm + 3.0_dp)*(amm + 5.0_dp)/((amm + 2.0_dp)*(amm + 1.0_dp)*ak*rc2)* &
                  (lam - am*(am + 1.0_dp) - ak*rc2/(amm + 3.0_dp))
         i1 = 4
      else
         rda(3) = (amm + 5.0_dp)*(amm + 7.0_dp)/((amm + 3.0_dp)*(amm + 2.0_dp)*ak*rc2)* &
                  (lam - (am + 1.0_dp)*(am + 2.0_dp) - (6.0_dp*am + 3.0_dp)*ak*rc2/((amm + 1.0_dp)* &
                  (amm + 5.0_dp)))
         i1 = 5
      end if
      do i = i1, nm2, 2
         ajj = i - 2.0_dp
         aii = 2.0_dp*ajj + amm
         if (guard) then
            den = rda(i-2) + 1.0e-50_dp
         else
            den = rda(i-2)
         end if
         rda(i) = (aii + 3.0_dp)*(aii + 5.0_dp)/((amm + ajj + 2.0_dp)*(amm + ajj + 1.0_dp))* &
                  ((lam - (am + ajj)*(am + ajj + 1.0_dp) - (2.0_dp*(am + ajj)*(am + ajj + 1.0_dp) &
                  - 2.0_dp*m2 - 1.0_dp)*ak*rc2/((aii - 1.0_dp)*(aii + 3.0_dp)))/ &
                  rc2*ak - ajj*(ajj - 1.0_dp)/((aii - 3.0_dp)*(aii - 1.0_dp))/den)
      end do
      raa = rda(nm2)

      ! iteration 'from top to below'
      rda(nmax) = czero
      do i = 2, nmax - 2, 2
         j = nmax - i
         aj = j
         ajj = amm + 2.0_dp*aj
         rda(j) = aj*(aj - 1.0_dp)*(ajj + 3.0_dp)*(ajj + 5.0_dp)/((ajj - 3.0_dp)*(ajj - 1.0_dp)* &
                  (amm + aj + 2.0_dp)*(amm + aj + 1.0_dp))/((ajj + 3.0_dp)*(ajj + 5.0_dp)/ &
                  ((amm + aj + 2.0_dp)*(amm + aj + 1.0_dp)*ak*rc2)*(lam - (am + aj)* &
                  (am + aj + 1.0_dp) - (2.0_dp*(am + aj)*(am + aj + 1.0_dp) - 2.0_dp*m2 - 1.0_dp)*ak*rc2/ &
                  ((ajj - 1.0_dp)*(ajj + 3.0_dp))) - rda(j+2))
      end do
      sl = raa - rda(nm2)
   end function cf_residual

   ! Small-c approximation of lambda (two terms; used by o_l, o_s, n and p)
   function lam_small(n, m, rc2) result(lam)
      integer, intent(in) :: n, m
      complex(dp), intent(in) :: rc2
      complex(dp) :: lam
      real(dp) :: am, an3, an4
      integer :: n1, n2

      n1 = n - m
      n2 = n + m
      an3 = 2.0_dp*n - 1.0_dp
      an4 = an3 + 4.0_dp
      am = 4.0_dp*m*m - 1.0_dp
      lam = n*(n + 1.0_dp) - ak*rc2/2.0_dp*(am/an3/an4 - 1.0_dp) + (rc2*rc2)/2.0_dp* &
            ((n1 - 1.0_dp)*n1*(n2 - 1.0_dp)*n2/((an3 + 2.0_dp)*an3**3*an4) - (n1 + 1.0_dp)* &
            (n1 + 2.0_dp)*(n2 + 1.0_dp)*(n2 + 2.0_dp)/((an4 + 2.0_dp)*an4**3*(an3 + 2.0_dp)))
   end function lam_small

   ! Large-c asymptotic approximation of oblate lambda (without the
   ! exponentially small term)
   function lam_oblate_large(n, m, c) result(lam)
      integer, intent(in) :: n, m
      complex(dp), intent(in) :: c
      complex(dp) :: lam, rc2, rc4
      real(dp) :: q, q2, q4, q6, d, rm2

      if (mod(n - m, 2) /= 0) then
         q = n
      else
         q = n + 1.0_dp
      end if
      rc2 = c*c
      rc4 = rc2*rc2
      rm2 = real(m*m, dp)
      q2 = q*q
      q4 = q2*q2
      q6 = q4*q2
      d = q2 + 1.0_dp - rm2
      lam = -rc2 + 2.0_dp*c*q - d/2.0_dp - q*d/8.0_dp/c &
            - (4.0_dp*q2*(d + 1.0_dp) + d**2)/64.0_dp/rc2 &
            - q*(33.0_dp*q4 + 114.0_dp*q2 + 37.0_dp - 2.0_dp*rm2*(23.0_dp*q2 + 25.0_dp) &
            + 13*rm2**2)/512.0_dp/rc2/c &
            - (63.0_dp*q6 + 340.0_dp*q4 + 239.0_dp*q2 + 14.0_dp &
            - 10.0_dp*rm2*(10.0_dp*q4 + 23.0_dp*q2 + 3.0_dp) + rm2*rm2*(39.0_dp*q2 - 18.0_dp) &
            - 2.0_dp*rm2**3)/1024.0_dp/rc4 &
            - q*(527.0_dp*q6 + 4139.0_dp*q4 + 5221.0_dp*q2 + 1009.0_dp &
            - rm2*(939.0_dp*q4 + 3750.0_dp*q2 + 1591.0_dp) + rm2*rm2*(465.0_dp*q2 + 635.0_dp) &
            - 53.0_dp*rm2**3)/8192.0_dp/rc4/c
   end function lam_oblate_large

!=============================================================================
!  CDLAMo_l -- oblate eigenvalues (large c)
!=============================================================================
   subroutine cdlamo_l(vl, m, nn, c, eps, ier)
      integer, intent(in) :: m, nn
      complex(dp), intent(inout) :: vl(nn)
      complex(dp), intent(in) :: c
      real(dp), intent(in) :: eps
      integer, intent(out) :: ier

      character(*), parameter :: fmt2000 = '(14x,"******* CDLAMo_l ****** NMAX=",i4)'
      character(*), parameter :: fmt3002 = '(1x,"L",2f7.3,2i5,2x,3(f10.3,f10.3))'
      complex(dp) :: rl(1550), sl(1550)
      complex(dp) :: dx, sc, sdelta, rc2, rc4, r0, rdel, r_cont
      real(dp) :: eps1, ac, cc, af, am2, cn, an, delta, cc1, cc2, cff, cbb
      real(dp) :: ch, cr, crn, cc0, ck, q, q2, q4, q6, as, rm2
      integer :: inull, in, inout, nt, iv, n1001, n, n1, nm2, nmax, jn, num
      integer :: n2, n3, n4

      inull = 0
      cn = 0.0_dp
      r0 = czero
991   continue
      eps1 = eps
      ac = real(c, dp)
      cc = abs(c)
      af = 0.0_dp
      if (abs(ac - cc) > 1.0e-10_dp) af = 0.1_dp
      sc = cmplx(0.1_dp, af, dp)
      ier = 0
      if (ac < 0.0_dp) then
         ier = 3
         return
      end if
      vl(1:nn) = czero
      rl(1:550) = czero
      sl(1:550) = czero
      sdelta = czero

      rm2 = real(m*m, dp)
      am2 = 4.0_dp*m*m - 1.0_dp
      rc2 = c*c
      rc4 = rc2*rc2

      do in = 1, nn
         inout = 4
         nt = 0
         iv = 0

         dx = cmplx(0.25_dp, 0.25_dp*af, dp)
         if (ac*m > 400.0_dp) dx = cmplx(0.35_dp + (ac*m - 400.0_dp)/1000.0_dp, 0.35_dp*af, dp)

         if (in > 6) then
            if (real(vl(in-1), dp) < -4*ac) then
               dx = cmplx(0.1_dp, 0.1_dp*af, dp)
               if (ac > 60.0_dp) dx = cmplx(0.25_dp, 0.25_dp*af, dp)
               if (ac > 70.0_dp) dx = cmplx(0.3_dp, 0.3_dp*af, dp)
            end if
            if (ac > 40.0_dp .and. real(vl(in-2), dp) < 0.0_dp .and. &
                real(vl(in-1), dp) > -4*ac) &
               dx = cmplx(0.35_dp + (ac - 40.0_dp)/200.0_dp, 0.35_dp*af, dp)
            if (ac > 70.0_dp .and. abs(real(vl(in-2), dp) + real(vl(in-1), dp)) < cn) &
               dx = cmplx(0.35_dp, 0.35_dp*af, dp)
         end if

         n1001 = 0
         if (real(sc, dp) < 0.0_dp .and. n1001 == 0) sc = -sc
         n = m + in - 1
         n1 = n - m
         nm2 = n - m + 2
         nmax = int(nm2 + 2*cc + 2)
1113     continue
         jn = 0
         if (nmax > 500) write (lu, fmt2000) nmax
         if (mod(nmax - nm2, 2) /= 0) nmax = nmax + 1
         if (ac < 5.50_dp) goto 401

         if (n < 7 .or. in <= 3) goto 401
         if (in < 7) goto 401
         if (real(vl(in-1), dp) < 0.0_dp .and. m >= 10) goto 401

         if (real(vl(in-1), dp) < 3.0_dp*ac .and. m >= 15) then
            if (real(sc, dp) < 0.0_dp .and. n1001 == 0) sc = -sc
            goto 401
         end if

         cn = ac
         an = n
         if (an < cn) cn = an

         if (abs(real(vl(in-1), dp)) < cn .and. real(vl(in-2), dp) < cn &
             .or. abs(real(vl(in-2) + vl(in-1), dp)) < ac + m) then
            if (real(vl(in-1), dp) < 0.0_dp) rl(1) = czero
            if (real(vl(in-1), dp) > 0.0_dp) rl(1) = 2.0_dp*vl(in-1)
            if (real(vl(in-1), dp) > cn) rl(1) = vl(in-1) + ac/10.0_dp
            r0 = rl(1)
            if (n1001 > 1 .and. m < 5) dx = dx/1.25_dp
            if (n1001 /= 0 .and. m >= 5) dx = 2.0_dp*dx
            if (n1001 == 0) then
               if (real(dx, dp) < ac/100.0_dp) dx = ac/100.0_dp
               if (real(sc, dp) > 0.0_dp .and. ac < 60.0_dp .and. m < 5) sc = -sc
            end if
            iv = 100
            goto 2
         end if

         if (real(vl(in-1), dp) < cn .and. abs(real(vl(in-2) + vl(in-1), dp)) < cn &
             .or. real(vl(in-1), dp) < cn .and. abs(real(vl(in-2) + vl(in-3), dp)) < cn) then
            nt = 2
            goto 401
         end if

         if (real(vl(in-1), dp) > cn .and. real(vl(in-1) + vl(in-2), dp) > ac + m &
             .or. real(vl(in-2), dp) > 0.0_dp) then

            if (real(dx, dp) > 0.1_dp .and. abs(vl(in-1) - vl(in-2)) > ac + 2*m + 5) &
               dx = cmplx(0.1_dp, 0.1_dp*af, dp)

            if (abs(vl(in-2) - vl(in-3)) < cn .or. abs(vl(in-1) - vl(in-2)) < cn) then
               if (mod(n - m + 1, 2) == 0) r0 = vl(in-1) + abs(vl(in-2) - vl(in-3))
               if (mod(n - m + 1, 2) /= 0) r0 = vl(in-2) + 2.0_dp*ac
               iv = 1111
               if (n1001 == 0 .or. n1001 == 1) then
                  dx = cmplx(0.25_dp + ac*m/10000.0_dp, 0.25_dp*af, dp)
                  if (ac > 60.0_dp) dx = cmplx(0.3_dp + ac*m/10000.0_dp, 0.3_dp*af, dp)
                  if (ac > 75.0_dp) dx = cmplx(0.4_dp + ac*m/10000.0_dp, 0.4_dp*af, dp)
               end if
               goto 223
            end if

            r0 = vl(in-1) + abs(vl(in-1) - vl(in-2))

            if (abs(vl(in-1) - vl(in-2)) < 10.0_dp - m/15.0_dp) then
               if (ac >= 40.0_dp) then
                  if (ac*m >= 600.0_dp .or. real(vl(in-1), dp) > cn + 2*m) &
                     r0 = r0 + vl(in-2)/2.0_dp
                  if (ac*m < 600.0_dp .or. real(vl(in-1), dp) < cn + 2*m) &
                     r0 = r0 + vl(in-2)
               end if
               if (ac < 40.0_dp .and. m <= 15) then
                  if (real(vl(in-1), dp) > ac + m) r0 = r0 + vl(in-2)/2.0_dp
                  if (real(vl(in-1), dp) <= ac + m) r0 = r0 + vl(in-2)
               end if
            end if

223         continue
            if (ac >= 30.0_dp .and. m > 3) r0 = r0 + m/5.0_dp
            if (ac >= 16.0_dp .and. ac < 29.0_dp) r0 = r0 + 2.086_dp
            if (ac >= 30.0_dp .and. ac < 40.0_dp) r0 = r0 + (ac - 30.0_dp)/2.0_dp
            if (ac >= 30.0_dp .and. ac < 40.0_dp .and. real(vl(in-2), dp) < cn) r0 = r0 - m
            if (ac >= 40.0_dp .and. ac < 50.0_dp) r0 = r0 + n/10.0_dp - (ac - 40.0_dp)/3.0_dp - m/2.5_dp
            if (ac >= 50.0_dp .and. ac <= 60.0_dp) r0 = r0 + (ac - 50.0_dp)/3.0_dp
            if (ac > 60.0_dp .and. ac <= 70.0_dp) r0 = r0 + (ac - 60.0_dp)/3.0_dp
            if (ac > 70.0_dp .and. ac <= 80.0_dp) r0 = r0 + (ac - 70.0_dp)/4.0_dp
            if (ac > 80.0_dp .and. ac <= 90.0_dp) r0 = r0 - (ac - 80.0_dp)/8.0_dp + 3.5_dp
            if (ac > 90.0_dp .and. ac <= 100.0_dp) r0 = r0 - (ac - 90.0_dp)/8.0_dp + 3.5_dp
            if (ac > 100.0_dp) r0 = r0 - (ac - 100.0_dp)/8.0_dp + 3.5_dp
            if (ac > 110.0_dp) r0 = r0 + (ac - 110.0_dp)/8.0_dp + 0.5_dp
            if (ac > 50.0_dp .and. m > 15) r0 = r0 - n/40.0_dp
            iv = 11
            goto 2
         end if

         nt = 0
         goto 401

400      continue
         jn = jn + 1
         rl(1) = r0 - jn*dx*sc
         if (jn == 1001) goto 1112
         goto 22

401      continue
         iv = 401
         rl(1) = czero
         r0 = rl(1)
         if (cc > 0.5_dp) then
            n3 = 2*n - 1
            n4 = n3 + 4
            if (cc > (n + 3.0_dp)) then
               ! initial approximation (large c)
               rl(1) = lam_oblate_large(n, m, c)
               if (in > 2) then
                  if (real(vl(in-1), dp) < 0.0_dp .and. real(vl(in-2), dp) < 0.0_dp &
                      .and. real(rl(1), dp) > 0.0_dp) then
                     if (inull == 0 .or. ac < 83.0_dp .or. real(rl(1), dp) >= ac/2.0_dp) then
                        rl(1) = czero
                     else
                        rl(1) = vl(in-1)/2.0_dp
                        sdelta = czero
                     end if
                  end if
               end if

               if (ac > 24.92_dp .and. ac < 24.93_dp .and. m == 1 .and. n == 15) then
                  rl(1) = -2.78_dp
                  inout = 32
               end if
               if (ac > 55.32_dp .and. ac < 55.33_dp .and. m == 1 .and. n == 34) then
                  rl(1) = -52.247_dp
                  inout = 32
                  sdelta = czero
               end if
               if (ac > 63.69_dp .and. ac < 63.70_dp .and. m == 1 .and. n == 37) then
                  rl(1) = -144.2_dp
                  inout = 32
                  sdelta = czero
               end if
               r0 = rl(1)
            else
               ! initial approximation (small c)
               rl(1) = lam_small(n, m, rc2)
               r0 = rl(1)
            end if
         end if

2        continue
         if (jn == 0) then
            delta = real(sdelta, dp)
            if (real(r0, dp) < 0.0_dp .or. nt /= 0) r0 = r0 + delta
            if (m > 6 .and. abs(real(r0, dp)) < cn .and. real(dx, dp) < 0.7_dp) dx = 0.7_dp
            if (m > 6 .and. n1001 == 0 .and. in > 1) then
               if (abs(real(r0, dp)) < cn .or. abs(real(vl(in-1), dp)) < cn) then
                  dx = 0.7_dp
                  if (m > 11 .and. ac*m > 400.0_dp .and. mod(n - m + 1, 2) /= 0) then
                     if (real(sc, dp) > 0.0_dp) sc = -sc
                     dx = 2.0_dp*(ac + m - 11)/10.0_dp
                  end if
               end if
            end if
            if (in > 1) then
               if (real(r0, dp) < real(vl(in-1), dp)) r0 = vl(in-1)
            end if
            rl(1) = r0
         end if

22       continue
         rdel = cmplx(0.1_dp, 0.1_dp*af, dp)
222      continue
         num = 1
         do
            sl(num) = cf_residual(rl(num), m, n1, nm2, nmax, rc2, .true.)
            num = num + 1
            if (num > 540) eps1 = 1.0e-10_dp
            if (mod(num, 2) /= 0) then
               rl(num) = rl(num-2)
               if (abs(sl(num-1) - sl(num-2)) <= eps1) exit
               rl(num) = rl(num-2) + sl(num-2)*rdel/(sl(num-2) - sl(num-1))
               if (abs(rl(num) - rl(num-2)) <= eps1) exit
            else
               if (num /= 2) rdel = rdel*sl(num-1)/(sl(num-3) - sl(num-2))
               rl(num) = rl(num-1) + rdel
            end if
         end do

         ! eigenvalue is found
         vl(in) = rl(num)
         if (ac < 2.5_dp) goto 110
         if (in <= 2) goto 110

         cc2 = real(rl(1), dp)
         cc1 = real(vl(in), dp)
         cff = abs(cc1 - cc2)
         if (cff > 20.0_dp*ac .and. cc1 < -ac) then
            rdel = rdel/2.0_dp
            rl(1) = r0
            goto 222
         end if
         cn = ac
         an = n
         if (an < cn) cn = an
         if (n1001 == 2 .and. real(vl(in), dp) > 0.0_dp .and. &
             real(vl(in-1), dp) < 0.0_dp .and. real(vl(in-2), dp) < 0.0_dp) cn = ac + an
         cbb = real(vl(in), dp) - real(vl(in-1), dp)

         if (cff > cn) goto 400
         if (cbb <= 0.0_dp .and. real(vl(in), dp) > -5*ac) goto 400
         if (cbb < 0.0_dp .and. real(vl(in), dp) <= -5*ac) goto 400

1112     continue
         ch = real(r0, dp)
         cr = real(rl(1), dp)
         crn = real(vl(in-2), dp)
         cc0 = real(vl(in-1), dp)
         ck = real(vl(in), dp)
         if (jn > 1) write (lu, fmt3002) c, n, jn, ch, cr, crn, cc0, ck
         if (jn > 1) write (*, fmt3002) c, n, jn, ch, cr, crn, cc0, ck
         if (jn < 1001) goto 110
         sc = -sc
         n1001 = n1001 + 1
         if (n1001 >= inout) then
            write (lu, *) ac
            write (*, *) ac
            write (lu, *) '**   Problems with computations of eigenvalues!'
            write (*, *) '**   Problems with computations of eigenvalues!'
            ier = 8
            return
         end if
         if (mod(n1001, 2) == 0 .and. nt == 0) dx = dx/3.0_dp
         goto 1113

110      continue
         if (real(vl(in), dp) > 0.0_dp) sdelta = vl(in) - r0
         if (real(vl(in), dp) < 0.0_dp .and. mod(in, 2) /= 0) sdelta = vl(in) - r0

         if (abs(sdelta) > 5.0_dp*ac .and. in < 3) then
            rdel = rdel/2.0_dp
            eps1 = eps1*10.0_dp
            write (*, *) in, rdel
            write (*, *) r0, vl(in)
            goto 222
         end if

         if (in > 2) then
            if (real(vl(in) - vl(in-1), dp) > in .or. &
                real(vl(in-1) - vl(in-2), dp) > in - 1 - ac/10.0_dp) cycle
         end if

         if (cc > n + 3.0_dp .and. iv /= 401) then
            if (mod(n1, 2) /= 0) then
               q = n
               as = 1.0_dp
            else
               q = n + 1.0_dp
               as = -1.0_dp
            end if
            q2 = q*q
            q4 = q2*q2
            q6 = q4*q2
            rl(1) = lam_oblate_large(n, m, c) &
                    + as*(4.0_dp*c)**(q + 1.0_dp)/exp(2.0_dp*c)/fact2(n + m, n) &
                    *(1.0_dp - (3.0_dp*q2 + 4.0_dp*q + 1.0_dp - rm2)/8.0_dp/c &
                    + ((9.0_dp*q4 + 4.0_dp*q2*q - 18.0_dp*q2 - 12.0_dp*q - 7.0_dp) &
                    - rm2*(6.0_dp*q2 - 4.0_dp*q - 6.0_dp) + rm2*rm2)/128.0_dp/rc2 &
                    - ((23.0_dp*q6 - 44.0_dp*q4*q - 171.0_dp*q4 + 400.0_dp*q2*q + 257.0_dp*q2 + 444.0_dp*q &
                    + 51.0_dp) - rm2*(19.0_dp*q4 - 56.0_dp*q2*q - 50.0_dp*q2 + 240.0_dp*q + as*55.0_dp) &
                    + rm2*rm2*(5.0_dp*q2 - 12.0_dp*q + 5.0_dp) - rm2*rm2*rm2)/3.0_dp/1024.0_dp/rc2/c)

            if (real(vl(in-1), dp) < 0.0_dp .and. real(vl(in-2), dp) < 0.0_dp &
                .and. real(rl(1), dp) > 0.0_dp) rl(1) = czero

            r0 = rl(1)
            if (abs(real(vl(in) - r0, dp)) > cn) then
               iv = 401
               if (real(sc, dp) < 0.0_dp) sc = -sc
               goto 2
            end if
         end if
      end do

      ! consistency check of the last eigenvalue against the small-c formula
      n = m + nn - 1
      n1 = n - m
      n3 = 2*n - 1
      n4 = n3 + 4
      n2 = n + m
      r_cont = n*(n + 1.0_dp) - ak*rc2/2.0_dp*(am2/n3/n4 - 1.0_dp) + rc4/2.0_dp* &
               ((n1 - 1.0_dp)*n1*(n2 - 1.0_dp)*n2/((n3 + 2.0_dp)*n3**3*n4) - (n1 + 1.0_dp)* &
               (n1 + 2.0_dp)*(n2 + 1.0_dp)*(n2 + 2.0_dp)/((n4 + 2.0_dp)*n4**3*(n3 + 2.0_dp)))

      if (abs(r_cont - vl(nn)) > 3.0_dp*ac) then
         inull = inull + 1
         write (lu, *) '**   ', r_cont
         write (lu, *) '**   ', vl(nn)
         write (*, *) '**   ', r_cont
         write (*, *) '**   ', vl(nn)
         if (inull == 2) then
            ier = 8
            return
         end if
         goto 991
      end if
   end subroutine cdlamo_l

!=============================================================================
!  CDLAMo_s -- oblate eigenvalues (small c)
!=============================================================================
   subroutine cdlamo_s(vl, m, nn, c, eps, ier)
      integer, intent(in) :: m, nn
      complex(dp), intent(inout) :: vl(nn)
      complex(dp), intent(in) :: c
      real(dp), intent(in) :: eps
      integer, intent(out) :: ier

      character(*), parameter :: fmt2000 = '(14x,"******* CDLAMo_s ****** NMAX=",i4)'
      character(*), parameter :: fmt3002 = '(1x,"s",2f7.3,2i5,2x,3(f10.3,f10.3))'
      complex(dp) :: rl(550), sl(550)
      complex(dp) :: dx, sc, rc2, r0, rdel
      real(dp) :: eps1, ac, cc, af, cn, an, cc1, cc2, cff, cbb, ch, cr, crn, cc0, ck
      integer :: in, inout, n1001, n, n1, nm2, nmax, jn, num, nc

      eps1 = eps
      ac = real(c, dp)
      cc = abs(c)
      af = 0.0_dp
      if (abs(ac - cc) > 1.0e-10_dp) af = 0.1_dp
      sc = cmplx(0.1_dp, af, dp)
      ier = 0
      if (ac < 0.0_dp) then
         ier = 3
         return
      end if
      vl(1:nn) = czero
      rl = czero
      sl = czero
      r0 = czero

      rc2 = c*c
      nc = int(ac)

      do in = 1, nn
         inout = 8
         dx = cmplx(0.1_dp, 0.1_dp*af, dp)
         if (ac > 20.0_dp .and. m > 10) dx = cmplx(0.1_dp + (m - 10.0_dp)/100.0_dp, 0.1_dp*af, dp)
         n1001 = 0
1113     continue
         jn = 0
         n = m + in - 1
         n1 = n - m
         nm2 = n - m + 2
         nmax = int(nm2 + 2*cc + 2)
         if (nmax > 500) write (lu, fmt2000) nmax
         if (mod(nmax - nm2, 2) /= 0) nmax = nmax + 1

         if (ac < 5.50_dp) goto 401
         if (in < 7) goto 401
         if (n < nc - nc/3) goto 401
         if (real(vl(in-1), dp) < 0.0_dp) goto 401
         ! (as in the original, rl(1) is not reset on this path)
         r0 = vl(in-1) + 2.0_dp*n + 1.0_dp
         if (ac > 18.0_dp) r0 = r0 - 2.532_dp

         if (m == 1 .and. ac > 15.998_dp .and. ac < 16.003_dp .and. in == 13) then
            r0 = 66.1923_dp - (ac - 16.0_dp)*10.0_dp
            dx = dx/100.0_dp
         end if
         if (m == 7 .and. ac > 18.690_dp .and. ac < 18.695_dp .and. in == 10) then
            r0 = 134.161_dp - (ac - 18.691_dp)*13.69744_dp
            dx = dx/5.5_dp
         end if
         if (m == 8 .and. ac > 18.535_dp .and. ac < 18.540_dp .and. in == 8) then
            r0 = 116.313_dp - (ac - 18.535_dp)*13.69744_dp
            dx = dx/5.5_dp
         end if
         goto 2

400      continue
         jn = jn + 1
         rl(1) = r0 - jn*dx*sc
         if (jn == 1001) goto 1112
         goto 2

401      continue
         rl(1) = czero
         if (cc > 0.5_dp) then
            if (cc > (n + 3.0_dp)) then
               rl(1) = lam_oblate_large(n, m, c)   ! large c
            else
               rl(1) = lam_small(n, m, rc2)        ! small c
            end if
            r0 = rl(1)
         end if

2        continue
         if (in > 1) then
            if (real(vl(in-1), dp) < 0.0_dp .and. jn == 1) then
               r0 = r0 - 2.0_dp
               if (real(r0, dp) > 0.0_dp) r0 = r0 - 2.0_dp
            end if
         end if
         rdel = dx
         num = 1
         do
            sl(num) = cf_residual(rl(num), m, n1, nm2, nmax, rc2, .false.)
            num = num + 1
            if (num > 540) eps1 = 1.0e-10_dp
            if (mod(num, 2) /= 0) then
               rl(num) = rl(num-2)
               if (abs(sl(num-1) - sl(num-2)) <= eps1) exit
               rl(num) = rl(num-2) + sl(num-2)*rdel/(sl(num-2) - sl(num-1))
               if (abs(rl(num) - rl(num-2)) <= eps1) exit
            else
               if (num /= 2) rdel = rdel*sl(num-1)/(sl(num-3) - sl(num-2))
               rl(num) = rl(num-1) + rdel
            end if
         end do

         ! eigenvalue is found
         vl(in) = rl(num)
         if (ac < 2.5_dp) cycle
         if (in <= 2) cycle
         cc2 = real(rl(1), dp)
         cc1 = real(vl(in), dp)
         cff = abs(cc1 - cc2)
         cn = ac
         an = n
         if (an < cn) cn = an
         if (an > cn .and. ac < 10.0_dp) cn = an
         cbb = real(vl(in), dp) - real(vl(in-1), dp)
         if (cff > cn .or. cbb <= 0.0_dp) goto 400

1112     continue
         ch = real(r0, dp)
         cr = real(rl(1), dp)
         crn = real(vl(in-2), dp)
         cc0 = real(vl(in-1), dp)
         ck = real(vl(in), dp)
         if (jn > 1) write (lu, fmt3002) c, n, jn, ch, cr, crn, cc0, ck
         if (jn > 1) write (*, fmt3002) c, n, jn, ch, cr, crn, cc0, ck
         if (jn < 1001) cycle

         sc = -sc
         n1001 = n1001 + 1
         if (n1001 >= inout) then
            write (lu, *) '**   Problems with computations of eigenvalues!'
            write (*, *) '**   Problems with computations of eigenvalues!'
            ier = 8
            return
         end if
         if (mod(n1001, 2) == 0) dx = dx/3.0_dp
         goto 1113
      end do
   end subroutine cdlamo_s

!=============================================================================
!  CDLAMp -- prolate eigenvalues
!=============================================================================
   subroutine cdlamp(vl, m, nn, c, eps, ier)
      integer, intent(in) :: m, nn
      complex(dp), intent(inout) :: vl(nn)
      complex(dp), intent(in) :: c
      real(dp), intent(in) :: eps
      integer, intent(out) :: ier

      character(*), parameter :: fmt2000 = '(14x,"******* CDLAMp ****** NMAX=",i4)'
      character(*), parameter :: fmt3002 = '(1x,"p",2f6.2,2i5,2x,3(f10.3,f10.3))'
      complex(dp) :: rl(1550), sl(1550)
      complex(dp) :: dx, sc, sdelta, rc2, rc4, r0, rdel
      real(dp) :: eps1, ac, cc, af, cn, an, delta, cc1, cc2, cff, cr, crn, cc0
      real(dp) :: rm2, rm4, q, q2, q3, q4, an3, an4
      integer :: in, inout, nt, n1001, n, n1, n2, nm2, nmax, jn, num, nc

      eps1 = eps
      ac = real(c, dp)
      cc = abs(c)
      af = 0.0_dp
      if (abs(ac - cc) > 1.0e-10_dp) af = 0.1_dp
      sc = cmplx(0.1_dp, af, dp)
      if (ac < 29.0_dp) sc = -sc
      ier = 0
      if (ac < 0.0_dp) then
         ier = 3
         return
      end if
      vl(1:nn) = czero
      rl = czero
      r0 = czero
      nt = 0
      sdelta = czero
      rm2 = real(m*m, dp)
      rm4 = rm2*rm2
      rc2 = c*c
      rc4 = rc2*rc2
      nc = int(ac)

      do in = 1, nn
         inout = 16
         dx = cmplx(0.1_dp, 0.1_dp*af, dp)
         if (ac > 60.0_dp) dx = cmplx(0.2_dp, 0.2_dp*af, dp)
         n1001 = 0
1113     continue
         jn = 0
         n = m + in - 1
         n1 = n - m
         nm2 = n - m + 2
         nmax = int(nm2 + 2*cc + 2)
         if (nmax > 3090) write (lu, fmt2000) nmax
         if (nmax > 3090) write (*, fmt2000) nmax
         if (mod(nmax - nm2, 2) /= 0) nmax = nmax + 1
         if (ac < 11.90_dp) goto 401
         if (in < 7) goto 401
         if (n < nc - nc/3 - m .and. nt == 0) goto 401
         ! (as in the original, rl(1) is not reset on this path)
         r0 = vl(in-1) + abs(real(vl(in-1), dp) - real(vl(in-2), dp))
         if (ac >= 30.0_dp .and. m > 3) r0 = r0 + m/5.0_dp
         if (ac >= 16.0_dp .and. ac < 29.0_dp) r0 = r0 + 1.5_dp
         if (ac >= 30.0_dp .and. ac < 40.0_dp) r0 = r0 + (ac - 30.0_dp)/2.0_dp
         if (ac >= 40.0_dp .and. ac < 50.0_dp) r0 = r0 + n/10.0_dp - (ac - 40.0_dp)/3.0_dp
         if (ac >= 50.0_dp .and. ac <= 60.0_dp) r0 = r0 + (ac - 50.0_dp)/3.0_dp
         if (ac > 60.0_dp .and. ac <= 70.0_dp) r0 = r0 + (ac - 60.0_dp)/3.0_dp
         if (ac > 70.0_dp .and. ac <= 80.0_dp) r0 = r0 + (ac - 70.0_dp)/4.0_dp
         if (ac > 80.0_dp .and. ac <= 90.0_dp) r0 = r0 - (ac - 80.0_dp)/8.0_dp + 3.5_dp
         if (ac > 90.0_dp .and. ac <= 100.0_dp) r0 = r0 - (ac - 90.0_dp)/8.0_dp + 4.0_dp
         if (ac > 100.0_dp) r0 = r0 - (ac - 100.0_dp)/8.0_dp + 4.5_dp
         if (ac > 110.0_dp) r0 = r0 + (ac - 110.0_dp)/8.0_dp + 1.5_dp
         goto 2

400      continue
         sdelta = czero
         jn = jn + 1
         rl(1) = r0 - jn*dx*sc*ak
         if (jn == 1001) goto 1112
         goto 2

401      continue
         rl(1) = czero
         r0 = czero
         if (cc > 4.3_dp) then
            if (cc > (n + 3.0_dp)) then
               ! initial approximation (large c)
               q = 2.0_dp*n1 + 1.0_dp
               q2 = q*q
               q3 = q2*q
               q4 = q3*q
               rl(1) = c*q + rm2 - (q2 + 5.0_dp)/8.0_dp - q*(q2 + 11.0_dp - 32.0_dp*rm2)/c/64.0_dp &
                       - (5.0_dp*(q2**2 + 26.0_dp*q2 + 21.0_dp) - 384.0_dp*rm2*(q2 + 1.0_dp))/1024.0_dp/rc2 &
                       - 1.0_dp/rc2/c*(1.0_dp/128.0_dp**2*(33.0_dp*q4*q + 1594.0_dp*q3 &
                       + 5621.0_dp*q) - rm2/128.0_dp*(37.0_dp*q3 + 167.0_dp*q) + rm4/8.0_dp*q) &
                       - 1.0_dp/rc2/rc2*(1.0_dp/256.0_dp**2*(63.0_dp*q4*q2 + 4940.0_dp*q4 &
                       + 43327.0_dp*q2 + 22470.0_dp) - rm2/512.0_dp*(115.0_dp*q4 &
                       + 1310.0_dp*q2 + 735.0_dp) + 3.0_dp*rm4/8.0_dp*(q2 + 1.0_dp)) &
                       - 1.0_dp/rc2/rc2/c*(1.0_dp/1024.0_dp**2*(527.0_dp*q4*q3 + 61529.0_dp*q4*q &
                       + 1043961.0_dp*q3 + 2241599.0_dp*q) - rm2/32.0_dp/1024.0_dp*(5739.0_dp*q4*q &
                       + 127550.0_dp*q3 + 298951.0_dp*q) + rm4/512.0_dp*(355.0_dp*q3 + 1505.0_dp*q) &
                       - rm4*rm2*q/16.0_dp)
            else
               ! initial approximation (small c): two terms + c^6 correction
               an3 = 2.0_dp*n - 1.0_dp
               an4 = an3 + 4.0_dp
               n2 = n + m
               rl(1) = lam_small(n, m, rc2) &
                       + rc4*rc2*(4.0_dp*rm2 - 1.0_dp)*((n1 + 1.0_dp)*(n1 + 2.0_dp)*(n2 + 1.0_dp)*(n2 + 2.0_dp) &
                       /(an3*(2.0_dp*n + 1.0_dp)*an4**5*(an4 + 2.0_dp)*(an4 + 4.0_dp)) &
                       - (n1 - 1.0_dp)*n1*(n2 - 1.0_dp)*n2 &
                       /((an3 - 4.0_dp)*(an3 - 2.0_dp)*an3**5*(an3 + 2.0_dp)*an4))
            end if
            r0 = rl(1)
         end if

2        continue
         rdel = 0.1_dp*sc
         delta = real(sdelta, dp)
         if (delta > 1.0_dp) delta = 1.0_dp
         if (delta < -1.0_dp) delta = -1.0_dp
         if (n1001 == 0) r0 = r0 + delta
         if (abs(sdelta) > 0.8_dp .and. nt == 0) then
            sdelta = czero
            nt = 1
         end if

         num = 1
         do
            sl(num) = cf_residual(rl(num), m, n1, nm2, nmax, rc2, .false.)
            num = num + 1
            if (num > 540) eps1 = 1.0e-10_dp
            if (mod(num, 2) /= 0) then
               if (num >= 1009) exit
               rl(num) = rl(num-2)
               if (abs(sl(num-1) - sl(num-2)) <= eps1) exit
               rl(num) = rl(num-2) + sl(num-2)*rdel/(sl(num-2) - sl(num-1))
               if (abs(rl(num) - rl(num-2)) <= eps1) exit
            else
               if (num /= 2) rdel = rdel*sl(num-1)/(sl(num-3) - sl(num-2))
               rl(num) = rl(num-1) + rdel
            end if
         end do

         ! eigenvalue is found
         vl(in) = rl(num)
         sdelta = vl(in) - r0
         if (ac < 11.90_dp) cycle
         if (in < 3) cycle
         cc2 = real(r0, dp)
         cc1 = real(vl(in), dp)
         cff = abs(cc1 - cc2)
         cn = ac
         an = n
         if (an < cn) cn = an
         if (cff >= cn) goto 400

1112     continue
         if (jn > 1) then
            cr = real(rl(1), dp)
            crn = real(vl(in-2), dp)
            cc0 = real(vl(in-1), dp)
            write (lu, fmt3002) c, n, jn, cc2, cr, crn, cc0, cc1
            write (*, fmt3002) c, n, jn, cc2, cr, crn, cc0, cc1
         end if
         if (jn < 1001) cycle
         sc = -sc
         n1001 = n1001 + 1
         if (mod(n1001, 2) == 0) dx = dx/1.5_dp

         if (n1001 >= inout) then
            write (lu, *) ac
            write (*, *) ac
            write (lu, *) '**   Problems with computations of eigenvalues!'
            write (*, *) '**   Problems with computations of eigenvalues!'
            ier = 8
            return
         end if
         goto 1113
      end do
   end subroutine cdlamp

!=============================================================================
!  CDLAMn -- eigenvalues at c, starting from the eigenvalues vl0 at a
!            nearby value of c (continuation step)
!=============================================================================
   subroutine cdlamn(vl0, vl, m, nn, c, eps, ier)
      integer, intent(in) :: m, nn
      complex(dp), intent(inout) :: vl0(nn), vl(nn)
      complex(dp), intent(in) :: c
      real(dp), intent(in) :: eps
      integer, intent(out) :: ier

      character(*), parameter :: fmt2000 = '(14x,"******* CDLAMn ****** NMAX=",i4)'
      character(*), parameter :: fmt3002 = '(1x,"*n*",2f6.2,2i5,2x,3(f10.3,f10.3),"*")'
      complex(dp) :: rl(1550), sl(1550)
      complex(dp) :: dx, sc, rc2, r0, rdel
      real(dp) :: eps1, ac, cc, af, cn, an, cc1, cc2, cff, cbb, ch, cr, crn, cc0, ck
      integer :: in, in0, n1001, n, n1, nm2, nmax, jn, num

      eps1 = eps
      ac = real(c, dp)
      cc = abs(c)
      af = 0.0_dp
      if (abs(ac - cc) > 1.0e-10_dp) af = 0.1_dp
      sc = cmplx(0.1_dp, af, dp)
      if (ac < 32 + m - 3 .or. ac > 47.4_dp) sc = -sc
      ier = 0
      in0 = 0
      if (ac < 0.0_dp) then
         ier = 3
         return
      end if
      vl(1:nn) = czero
      rc2 = c*c

      do in = 1, nn
         dx = cmplx(0.1_dp, 0.1_dp*af, dp)
         n1001 = 0
1113     continue
         jn = 0
         n = m + in - 1
         n1 = n - m
         nm2 = n - m + 2
         nmax = int(nm2 + 2*cc + 2)
         if (nmax > 1550) write (lu, fmt2000) nmax
         if (nmax > 1550) write (*, fmt2000) nmax
         if (mod(nmax - nm2, 2) /= 0) nmax = nmax + 1

         if (abs(vl0(in)) < 1.0e-10_dp) then
            in0 = in
            ! (guarded: the original read vl0(0), vl0(-1) when in <= 2)
            if (in > 2) vl0(in) = vl0(in-1) + abs(real(vl0(in-1) - vl0(in-2), dp))
            if (cc <= (n + 3.0_dp)) vl0(in) = lam_small(n, m, rc2)   ! small c
         end if

         r0 = vl0(in)
         rl(1) = r0
         goto 2

400      continue
         jn = jn + 1
         rl(1) = r0 - jn*dx*sc*ak
         if (jn == 1001 .and. in > 2) goto 1112
         if (jn == 1001 .and. in <= 2) goto 1131

2        continue
         rdel = dx
         num = 1
         do
            sl(num) = cf_residual(rl(num), m, n1, nm2, nmax, rc2, .false.)
            num = num + 1
            if (num > 1540) eps1 = eps1*10.0_dp
            if (num == 1550) then
               ier = 8
               return
            end if
            if (mod(num, 2) /= 0) then
               if (num >= 1009) exit
               rl(num) = rl(num-2)
               if (abs(sl(num-1) - sl(num-2)) <= eps1) exit
               rl(num) = rl(num-2) + sl(num-2)*rdel/(sl(num-2) - sl(num-1))
               if (abs(rl(num) - rl(num-2)) <= eps1) exit
            else
               if (num /= 2) rdel = rdel*sl(num-1)/(sl(num-3) - sl(num-2))
               rl(num) = rl(num-1) + rdel
            end if
         end do

         ! eigenvalue is found
         vl(in) = rl(num)
         if (ac < 2.5_dp) cycle
         cc2 = real(r0, dp)
         cc1 = real(vl(in), dp)
         cff = abs(cc1 - cc2)

         cn = ac
         an = n
         if (ac > 13.0_dp) an = an + 1
         if (an < cn) cn = an
         if (an > cn .and. ac < 10.0_dp) cn = an

         if (in <= 2) goto 1121
         if (in0 /= 0) cycle

         cbb = real(vl(in), dp) - real(vl(in-1), dp)
         if (kob == 0 .or. kob == 1 .and. aimag(c) > 0.5_dp) then
            if (cbb < 0.0_dp) cycle
         end if
         if (cff > cn .or. cbb <= 0.0_dp) goto 400

1112     continue
         ch = real(r0, dp)
         cr = real(rl(1), dp)
         crn = real(vl(in-2), dp)
         cc0 = real(vl(in-1), dp)
         ck = real(vl(in), dp)
         if (jn > 1) write (lu, fmt3002) c, n, jn, ch, cr, crn, cc0, ck
         if (jn > 1) write (*, fmt3002) c, n, jn, ch, cr, crn, cc0, ck
         if (jn < 1001) cycle
         sc = -sc
         n1001 = n1001 + 1
         eps1 = eps1*10.0_dp
         if (n1001 >= 4) then
            ier = 8
            return
         end if
         goto 119      ! (n1001 <= 3 always holds here)

1121     continue
         if (cff > ac) goto 400
1131     continue
         ch = real(r0, dp)
         cr = real(rl(1), dp)
         ck = real(vl(in), dp)
         if (jn > 1) write (lu, fmt3002) c, n, jn, ch, cr, ck
         if (jn > 1) write (*, fmt3002) c, n, jn, ch, cr, ck
         if (jn < 1001) cycle
         sc = -sc
         n1001 = n1001 + 1
         eps1 = eps1*10.0_dp
         if (n1001 >= 4) then
            ier = 8
            return
         end if

119      continue
         if (mod(n1001, 2) == 0) dx = dx/1.5_dp
         goto 1113
      end do
   end subroutine cdlamn

!=============================================================================
!  Factorial helpers
!=============================================================================
   pure function fact2(n, k) result(f)      ! n! / k!
      integer, intent(in) :: n, k
      real(dp) :: f
      integer :: i
      f = 1.0_dp
      do i = k + 1, n
         f = f*i
      end do
   end function fact2

   pure function fact3(n, k, f1) result(f)  ! f1 * n! / k!
      integer, intent(in) :: n, k
      real(dp), intent(in) :: f1
      real(dp) :: f
      integer :: i
      f = f1
      do i = k + 1, n
         f = f*i
      end do
   end function fact3

!=============================================================================
!  Associated Legendre functions
!=============================================================================

   ! P_{m+l-1}^m and derivative-related values tau for |x| arbitrary
   subroutine dugolm1(x, m, lmax, pim, tau)
      real(dp), intent(in) :: x
      integer, intent(in) :: m, lmax
      real(dp), intent(out) :: pim(lmax), tau(lmax)
      real(dp) :: a
      integer :: m2, i, l1

      m2 = 2*m
      a = 1.0_dp
      do i = 3, m2 - 1
         if (mod(i, 2) == 1) a = a*i
      end do
      pim(1) = a
      pim(2) = a*x*(m2 + 1.0_dp)
      tau(1) = a*m*x
      tau(2) = a*(m2 + 1.0_dp)*((m + 1.0_dp)*x*x - 1.0_dp)
      do l1 = 3, lmax
         pim(l1) = ((2.0_dp*l1 + m2 - 3.0_dp)*x*pim(l1-1) - (l1 + m2 - 2.0_dp)*pim(l1-2)) &
                   /(l1 - 1.0_dp)
         tau(l1) = (l1 + m - 1.0_dp)*x*pim(l1) - (l1 + m2 - 1.0_dp)*pim(l1-1)
      end do
   end subroutine dugolm1

   ! Fills p_leg, pd_leg, q_leg, qd_leg (former COMMON /LEG/)
   subroutine funlegnn(m, x, nf)
      integer, intent(in) :: m, nf
      real(dp), intent(in) :: x
      real(dp) :: a1, d, dtmp
      integer :: i, ier1, ier2, ier3, ier4

      if (nf > 2*nterms) then
         write (lu, '(1x,"funlegnn  NF > 2*nterms, NF=",2i10)') nf
         return
      end if
      p_leg = 0.0_dp
      pd_leg = 0.0_dp
      q_leg = 0.0_dp
      qd_leg = 0.0_dp

      if (kob == 0) then
         call dugolm1(x, m, nf, p_leg, pd_leg)
         if (x < 1.0_dp) then
            a1 = sqrt(1.0_dp - x*x)
         else
            a1 = sqrt(x*x - 1.0_dp)
         end if
         do i = 1, nf
            p_leg(i) = p_leg(i)*a1**m
            pd_leg(i) = -pd_leg(i)*a1**(m - 2)
            if (x > 1.0_dp) pd_leg(i) = -pd_leg(i)
         end do
      end if

      ier1 = 0
      dtmp = 0.0_dp
      if (kob == 1) call dlegf1(p_leg(1), dtmp, x, 0.0_dp, 0.0_dp, m, m, kob, 0, ier1)
      dtmp = 0.0_dp
      call dlegf2(q_leg(1), dtmp, x, 0.0_dp, -m, m, kob, 0, ier2)
      if ((ier1 + ier2) /= 0) write (lu, '(1x,"funlegnn  M,I,IER1,2,3,4=",6i10)') m, ier1, ier2

      d = 0.0_dp
      ier3 = 0
      do i = 2, nf
         if (kob == 1) then
            if (i > 2) d = p_leg(i-2)
            call dlegf1(p_leg(i), pd_leg(i-1), x, p_leg(i-1), d, m + i - 1, m, kob, 1, ier3)
         end if
         call dlegf2(q_leg(i), qd_leg(i-1), x, q_leg(i-1), -m + i - 1, m, kob, 1, ier4)
         if (ier4 /= 0) write (lu, '(1x,"funlegnn  M,I,IER1,2,3,4=",6i10)') m, i, ier4
      end do
   end subroutine funlegnn

   ! P_n^m(0) (and derivative); not called by rad_fun, kept from the original
   subroutine dlegf00(res, red, n, m, id, ier)
      real(dp), intent(out) :: res, red
      integer, intent(in) :: n, m, id
      integer, intent(out) :: ier
      real(dp) :: e
      integer :: k1, k2, kk

      ier = 0
      res = 0.0_dp
      red = 0.0_dp
      if (n == 0 .and. m == 0) then
         res = 1.0_dp
         return
      end if
      if (n <= 0 .or. m < 0) ier = 1
      if (ier /= 0) return
      if (mod(abs(n - m), 2) == 0) then
         k1 = (m + n)/2
         e = 1.0_dp
         do kk = 1, k1
            e = e*(k1 + kk)
         end do
         k2 = (n - m)/2
         if (k2 + 1 > 170) error stop 1703
         res = (-1)**k2*e/fact(k2+1)/2.0_dp**n
      end if
      if (id == 0) return
      if (mod(abs(n - m), 2) /= 0) then
         k1 = (m + n + 1)/2
         e = 1.0_dp
         do kk = 1, k1
            e = e*(k1 + kk)
         end do
         k2 = (n - m - 1)/2
         red = (-1)**k2*e/fact(k2+1)/2.0_dp**n
      end if
   end subroutine dlegf00

   ! Associated Legendre function of the first kind (ic = 0: real argument,
   ! ic /= 0: imaginary argument), with derivative for id /= 0
   subroutine dlegf1(rs, rd, x1, x2, x3, n, m, ic, id, ier)
      real(dp), intent(inout) :: rs, rd
      real(dp), intent(in) :: x1, x2, x3
      integer, intent(in) :: n, m, ic, id
      integer, intent(out) :: ier

      real(dp) :: res, red, eps4, xx2, e, z, d, h, dd
      integer :: nm, j, jm, k1, k, k3, k4, i, isg

      ier = 0
      res = 0.0_dp
      red = 0.0_dp
      eps4 = 1.0e-20_dp
      xx2 = abs(x2)
      nm = n - m
      j = mod(n, 4)
      jm = mod(m, 4)
      if (n == 0 .and. m == 0) res = 1.0_dp
      if (n == 0 .and. m == 0) goto 100
      if (n <= 0 .or. m < 0) ier = 1
      if (id /= 0 .and. abs(x1 - 1.0_dp) < 1.0e-6_dp .and. ic == 0) ier = 2
      if (ier /= 0) goto 100
      if (abs(x1) > 1.0e-6_dp) goto 20

      ! ---- x1 = 0 ----
      if (mod(abs(nm), 2) /= 0) goto 100
      k1 = (m + n)/2
      e = 1.0_dp
      do k = 1, k1
         e = e*(k1 + k)
      end do
      if (nm/2 + 1 > 170) error stop 1703
      res = e/fact(nm/2+1)/2.0_dp**n
      if (id /= 0) goto 50
      goto 102

      ! ---- x1 /= 0 ----
20    continue
      if (m > n) goto 100
      if (nm >= 2 .and. x1 > 1.0_dp .and. ic == 0) goto 48
      z = 1.0_dp/(2.0_dp**n)
      k1 = (nm + 2)/2
      d = 0.0_dp
      if (ic == 0) then
         if (nm > 169) error stop 170
         do k = 1, k1
            k3 = nm - 2*(k - 1)
            if (k3 + 1 > 170) error stop 1704
            h = 1.0_dp/fact(k3+1)
            k4 = n - k + 1
            do i = 1, k4
               if (abs(h) > 1.0e300_dp) print *, n, h
               h = h*(k4 + i)
            end do
            ! (the original stored this sign in the argument ID)
            isg = 1
            if (mod(k - 1, 2) /= 0) isg = -1
            dd = isg*h/fact(k)*x1**k3
            d = d + dd
            if (abs(dd/d) < eps4) exit
         end do
         res = d*z*sqrt(abs(1.0_dp - x1**2)**m)
      else
         if (nm > 169) error stop 1705
         do k = 1, k1
            k3 = nm - 2*(k - 1)
            h = 1.0_dp/fact(k3+1)
            k4 = n - k + 1
            do i = 1, k4
               h = h*(k4 + i)
            end do
            d = d + h/fact(k)*x1**k3
         end do
         res = d*z*sqrt((1.0_dp + x1**2)**m)
      end if
      if (id /= 0) goto 50
      goto 102

48    continue
      res = ((2.0_dp*n - 1.0_dp)*x1*x2 - (n + m - 1.0_dp)*x3)/nm
      if (id /= 0) goto 50
      goto 102

      ! ---- derivative ----
50    continue
      if (ic == 0) then
         red = (nm*res - x1*n*x2)/(x1**2 - 1.0_dp)
         rs = res
         rd = red
         return
      end if
      red = (nm*res - x1*n*xx2)/(x1**2 - ak)

100   continue
      if (jm == 2 .or. jm == 3) then
         rs = -abs(res)
         if (j == 0 .or. j == 1) rs = -rs
      else
         rs = res
         if (j == 2 .or. j == 3) rs = -abs(res)
      end if
      j = mod(n - 1, 4)
      rd = red
      if (j == 2 .or. j == 3) rd = -abs(red)
      return

102   continue
      rs = res
      if (ic == 0 .and. jm == 3 .and. x1 > 1.0_dp) rs = abs(res)
      if (ic /= 0 .and. jm == 2) rs = -abs(res)
      if (ic /= 0 .and. jm == 3) rs = -abs(res)
   end subroutine dlegf1

   ! Associated Legendre function of the second kind, with derivative for
   ! id /= 0
   subroutine dlegf2(res, red, x1, x2, n, m, ic, id, ier)
      real(dp), intent(out) :: res
      real(dp), intent(inout) :: red
      real(dp), intent(in) :: x1, x2
      integer, intent(in) :: n, m, ic, id
      integer, intent(out) :: ier

      real(dp) :: r(1500), q(1500), d, e, q0, xer
      integer :: k1, k2, i, j, nmax

      ier = 0
      res = 0.0_dp
      red = 0.0_dp
      if (n + m < 0) ier = 1
      if (ic == 0 .and. abs(x1) <= 1.0_dp) ier = 2
      if (ic /= 0 .and. abs(x1) < 0.0_dp) ier = 2
      if (n + m - 1 < 0 .and. id /= 0) ier = 3
      if (ier /= 0) return
      if (abs(x1) > 1.0e-6_dp .or. ic == 0) goto 20

      ! ---- x1 = 0 and ic /= 0 ----
      if (m <= n) then
         k1 = n - m
         k2 = n + m - 1
         d = 1.0_dp
         e = 1.0_dp
         if (mod(k1, 2) == 0) then
            if (k1 /= 0) then
               do i = 1, k2, 2
                  d = d*i
               end do
               do i = 2, k1, 2
                  e = e*i
               end do
            end if
            res = -pi/2.0_dp*d/e
         else
            do i = 2, k2, 2
               d = d*i
            end do
            do i = 1, k1, 2
               e = e*i
            end do
            res = d/e
         end if
      else
         if (mod(abs(n - m), 2) == 0) return
         k1 = (n - m - 1)/2
         d = 1.0_dp
         do i = 1, k1
            d = d*(k1 + i)
         end do
         if (k1 + 1 > 170) error stop 1706
         res = 2.0_dp**n*d*fact(k1+1)
      end if
      goto 50

      ! ---- x1 /= 0 ----
20    continue
      if (m == 0) goto 61
      if (m == 1) goto 71

      ! m > 1
      if (ic == 0) res = log((x1 + 1.0_dp)/(x1 - 1.0_dp))/2.0_dp
      if (ic /= 0 .and. x1 < 0.5_dp) res = atan(x1) - pi/2.0_dp
      if (ic /= 0 .and. x1 >= 0.5_dp) res = -atan(1.0_dp/x1)
      r(1) = -ak/(sqrt(x1**2 - ak))
      r(2) = -2.0_dp*x1*ak/(sqrt(x1**2 - ak))*r(1)
      do i = 3, m
         r(i) = -2.0_dp*(i - 1.0_dp)*x1*ak/(sqrt(x1**2 - ak))*r(i-1) + &
                (2.0_dp - i)*(i - 1.0_dp)*r(i-2)
      end do
      q0 = r(m)
      if (n == 0) then
         res = q0
         goto 50
      end if
      j = n
      if (n < 0) j = -n
      nmax = int(30*x1 + j + 50)
      if (nmax > 500) nmax = 500
      r(nmax) = x1 - sqrt(x1**2 - ak)
      xer = abs(x1 - 1.0_dp/sqrt(3.0_dp))
      do i = nmax - 1, 1, -1
         r(i) = ((i + m)*ak)/((2.0_dp*i + 1.0_dp)*x1 - (i - m + 1.0_dp)*r(i+1))
      end do
      q(1) = r(1)*q0*ak
      if (xer < 1.0e-10_dp .and. m == 3) q(1) = 3.0_dp
      if (xer < 1.0e-10_dp .and. m == 6) q(1) = -144.0_dp
      if (xer < 1.0e-10_dp .and. m == 9) q(1) = 4.536e4_dp
      if (xer < 1.0e-10_dp .and. m == 12) q(1) = -4.35456e+07_dp
      if (xer < 1.0e-10_dp .and. m == 15) q(1) = 9.340529643126767e+10_dp
      if (xer < 1.0e-10_dp .and. m == 18) q(1) = -3.76610155207262e+14_dp
      if (xer < 1.0e-10_dp .and. m == 21) q(1) = 2.554546682742534e+18_dp
      if (xer < 1.0e-10_dp .and. mod(m, 3) == 0) q0 = 0.0_dp
      if (n == 1) then
         res = q(n)
         goto 50
      end if
      if (n > 0) then
         do i = 1, j - 1
            q(i+1) = r(i+1)*q(i)*ak**(i + 1)
         end do
         res = q(n)
         goto 50
      end if
      ! n < 0
      if (n == -1) then
         res = (ak*x1*q0 + (m - 1.0_dp)*q(1))/m
         goto 50
      end if
      r(1) = q0
      r(2) = (ak*x1*q0 + (m - 1.0_dp)*q(1))/m
      do i = -2, n, -1
         r(1-i) = (ak**i*(2.0_dp*i + 3.0_dp)*x1*r(-i) - &
                  (i - m + 2.0_dp)*r(-i-1))/(i + m + 1.0_dp)
      end do
      res = r(-n+1)
      goto 50

      ! m = 0
61    continue
      if (ic == 0) res = log((x1 + 1.0_dp)/(x1 - 1.0_dp))/2.0_dp
      if (ic /= 0 .and. x1 < 0.5_dp) res = atan(x1) - pi/2.0_dp
      if (ic /= 0 .and. x1 >= 0.5_dp) res = -atan(1.0_dp/x1)
      if (n /= 0) then
         q0 = res
         nmax = int(1.5_dp*x1 + n + 10)
         r(nmax) = x1 - sqrt(x1**2 - ak)
         do i = nmax - 1, 1, -1
            r(i) = ((i + m)*ak)/((2.0_dp*i + 1.0_dp)*x1 - (i - m + 1.0_dp)*r(i+1))
         end do
         q(1) = r(1)*q0*ak
         do i = 1, n - 1
            q(i+1) = r(i+1)*q(i)*ak**(i + 1)
         end do
         res = q(n)
      end if
      goto 50

      ! m = 1
71    continue
      res = -ak/(sqrt(x1**2 - ak))
      if (n /= 0) then
         q0 = res
         if (n == -1) then
            res = ak*x1*q0
         else
            nmax = int(30*x1 + n + 50)
            r(nmax) = x1 - sqrt(x1**2 - ak)
            do i = nmax - 1, 1, -1
               r(i) = ((i + m)*ak)/((2.0_dp*i + 1.0_dp)*x1 - (i - m + 1.0_dp)*r(i+1))
            end do
            q(1) = r(1)*q0*ak
            do i = 1, n - 1
               q(i+1) = r(i+1)*q(i)*ak**(i + 1)
            end do
            res = q(n)
         end if
      end if

      ! ---- derivative ----
50    continue
      if (id == 0) return
      if (ic == 0 .or. mod(n, 2) == 0) then
         red = ((n - m)*res - x1*n*x2)/(x1**2 - ak)
      else
         red = -((n - m)*res + x1*n*x2)/(x1**2 - ak)
      end if
   end subroutine dlegf2

end module spheroidal_radial

!=============================================================================
!  NOTES ON THE TRANSLATION
!
!  * COMMON /K1/, /EPS3/, /EPS33/, /FACT/, /PI/, /LEG/, /F2/ became module
!    variables; S and S1 became the constants czero and ci; K became kob.
!  * The Lahey-specific DREAL function was removed (intrinsic REAL(z,dp)).
!  * PAUSE (deleted from the standard) was replaced by a warning; the
!    computation continues as it did after the PAUSE.
!  * STOP nnnn  ->  ERROR STOP nnnn.
!  * Work arrays for Bessel functions and BD were enlarged by one element
!    (the original wrote one element past their end in the extreme case
!    num = 4*nterms).
!  * cdcof4a: an undefined variable (J3) in an overflow-rescaling loop was
!    replaced by L4, which is what the neighbouring loop uses.
!  * dlegf1: a sign was stored in the dummy argument ID (which may be a
!    constant); a local variable is used now.  This branch (ic = 0) is not
!    reached from rad_fun.
!  * cdrf2k02: the phase factor ss is now computed for all i12 (it was
!    undefined for i12 = 0, which rad_fun never uses).
!  * The RKF45 integrator keeps its state between calls; this is made
!    explicit with SAVE.
!  * Not changed, but worth checking (possible bugs in the original):
!      - rad_fun: W and W1 are REAL in the original (implicit typing), so
!        only the real part of the Wronskian is compared.
!      - cdrfp2: AKK = (Y1*Y3 + Y2*Y4*(1-G)/SUM) -- the parenthesisation
!        looks as if (Y1*Y3 + Y2*Y4)*(1-G)/SUM might have been intended.
!      - cdcof4a: after overflow rescaling the factor with index i is
!        applied twice.
!      - cdrf12cc: indices into the Legendre arrays can exceed 2*nterms
!        for large NK if the series does not converge early.
!=============================================================================
