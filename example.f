c     Example of using SPHFUN (see the header of sphfun.f)
c     compile:  gfortran -std=legacy -O2 example.f sphfun.f -o example
      PROGRAM EXAMPL
      PARAMETER (NMAX=330, LDD=330, LDS=19)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 C, RLAM(NMAX), DC(LDD,NMAX), RNORM(NMAX),
     *           R1(NMAX), R1D(NMAX), R2(NMAX), R2D(NMAX),
     *           SA(LDS,NMAX), SAD(LDS,NMAX)
      DIMENSION ETA(2), WERR(NMAX), METH(NMAX)
      DATA ETA /0.3D0, -0.7D0/
c
      CALL SPHINI
c
c.... prolate functions, m = 1, n = 1..10, real c = 5, ksi = 1.25:
c     radial functions of the 1st and 2nd kind and angular functions
      M  = 1
      NE = 10
      C  = (5D0, 0D0)
      XI = 1.25D0
      CALL SPHFUN(0,M,NE,C,XI,2,2,ETA,1D-10,0,
     *            RLAM,DC,LDD,RNORM,R1,R1D,R2,R2D,
     *            SA,SAD,LDS,WERR,METH,IER)
      WRITE(*,'(A,I3)') ' prolate, real c: IER =', IER
      WRITE(*,'(A)') '   n meth  lambda          R1              R1''
     *             R2              R2''             S(eta1)
     *    Wronsk.err'
      DO 10 L = 1, NE
      WRITE(*,'(2I4,6(1PE16.8),1PE10.2)') M+L-1, METH(L),
     *      DREAL(RLAM(L)), DREAL(R1(L)), DREAL(R1D(L)),
     *      DREAL(R2(L)), DREAL(R2D(L)), DREAL(SA(1,L)), WERR(L)
   10 CONTINUE
c
c.... same, but complex c (e.g. inside a particle with refractive
c     index 1.5+0.02i): 1st kind only
      C = (5D0, 0D0) * (1.5D0, 0.02D0)
      CALL SPHFUN(0,M,NE,C,XI,1,2,ETA,1D-10,0,
     *            RLAM,DC,LDD,RNORM,R1,R1D,R2,R2D,
     *            SA,SAD,LDS,WERR,METH,IER)
      WRITE(*,'(/A,I3)') ' prolate, complex c: IER =', IER
      DO 20 L = 1, NE
      WRITE(*,'(2I4,6(1PE16.8))') M+L-1, METH(L), RLAM(L), R1(L),
     *      R1D(L)
   20 CONTINUE
      END
