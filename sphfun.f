c=======================================================================
c
c     SPHFUN -- spheroidal wave functions (prolate and oblate, real or
c               complex size parameter c), extracted from program
c               hom6_6q ("Light Scattering by Spheroidal Particles",
c               N.V. Voshchinnikov and V.G. Farafonov, Astronomical
c               Institute, St. Petersburg University, (c) 2001-14).
c
c     The numerical routines below the driver are copied unchanged
c     from hom6_6q.  Only the driver SPHFUN (which reproduces the
c     logic of the subroutine HOMFUNQ of hom6_6q), the initialisation
c     routine SPHINI and the small routine CASWF (angular functions
c     with complex expansion coefficients, cf. DASWF1) are new.
c
c     Method: N.V. Voshchinnikov, V.G. Farafonov, Ap. Space Sci. 204,
c     19-86, 1993; Comput. Mathem. Mathem. Physics 43, 1299-1309, 2003.
c
c     Notation (as in hom6_6q):
c       prolate (KTYPE=0):  radial coordinate ksi > 1,
c         d/dksi[(ksi^2-1) dR/dksi] - (lambda - c^2 ksi^2
c                                      + m^2/(ksi^2-1)) R = 0,
c         d/deta[(1-eta^2) dS/deta] + (lambda - c^2 eta^2
c                                      - m^2/(1-eta^2)) S = 0;
c       oblate  (KTYPE=1):  radial coordinate ksi >= 0,
c         d/dksi[(ksi^2+1) dR/dksi] - (lambda - c^2 ksi^2
c                                      - m^2/(ksi^2+1)) R = 0,
c         d/deta[(1-eta^2) dS/deta] + (lambda + c^2 eta^2
c                                      - m^2/(1-eta^2)) S = 0.
c
c       S_mn(c,eta) = sum'_r d_r^mn(c) P_{m+r}^m(eta),
c       r = 0,2,4,... (n-m even) or r = 1,3,5,... (n-m odd),
c       P_n^m(eta) = (1-eta^2)^(m/2) d^m P_n(eta)/deta^m
c       (no Condon-Shortley phase (-1)^m).
c       N_mn = 2 sum'_r (d_r^mn)^2 (r+2m)! / ((2r+2m+1) r!)
c       (no complex conjugation: N_mn is complex for complex c).
c       Radial functions of the 1st and 2nd kind have the asymptotic
c       normalisation R ~ cos/sin(c*ksi - (n+1)pi/2) / (c*ksi), so that
c       for real c the Wronskian is
c         R1*R2' - R1'*R2 = 1/(c (ksi^2-1))  (prolate),
c         R1*R2' - R1'*R2 = 1/(c (ksi^2+1))  (oblate).
c
c     Usage:
c       CALL SPHINI                    (once, before the first SPHFUN)
c       CALL SPHFUN(KTYPE,M,NE,C,KSI,IKIND,NETA,ETA,EPSW,NKIN,
c    *              RLAM,DCOEF,LDD,RNORM,R1,R1D,R2,R2D,
c    *              SA,SAD,LDS,WERR,METH,IER)
c
c     Input:
c       KTYPE  = 0 prolate, = 1 oblate spheroidal functions
c       M      azimuthal index, 1 <= M <= ~ 100 (hom6_6q uses M >= 1;
c              M = 0 is NOT supported)
c       NE     number of functions: n = M, M+1, ..., M+NE-1
c              (1 <= NE <= 330)
c       C      complex*16 size parameter, Re C > 0.  For a homogeneous
c              particle: C = c1 = k*f outside, C = c2 = m_r*c1 inside
c       KSI    real*8 radial coordinate (prolate: KSI > 1,
c              oblate: KSI >= 0)
c       IKIND  = 1 radial functions of the 1st kind only
c                  (any C; this is the branch hom6_6q uses for c2),
c              = 2 radial functions of the 1st and 2nd kind
c                  (only for real C, i.e. Im C = 0; this is the branch
c                  hom6_6q uses for c1)
c       NETA   number of angular arguments (0 <= NETA <= 19)
c       ETA    real*8 ETA(NETA), angular arguments, |ETA| < 1
c              (for |ETA| = 1, S and dS/deta are returned as 0)
c       EPSW   relative accuracy required for the Wronskian check
c              (IKIND = 2), e.g. 1d-10 (default if EPSW <= 0)
c       NKIN   number of expansion coefficients d_r to use;
c              NKIN <= 0 -> default as in hom6_6q:
c              max(NE+40,60) (+40 if KSI > 1.5), at most 329
c       LDD    leading dimension of DCOEF (LDD >= NK used, i.e.
c              LDD >= 329 is always safe)
c       LDS    leading dimension of SA, SAD (LDS >= max(NETA,1))
c
c     Output (l = 1..NE corresponds to n = M+l-1):
c       RLAM(l)       complex*16 eigenvalue lambda_mn(c)
c       DCOEF(i,l)    complex*16 coefficients d_r^mn, i = 1..NK,
c                     r = 2(i-1) (n-m even) or r = 2(i-1)+1 (n-m odd)
c       RNORM(l)      complex*16 norm N_mn (see above)
c       R1(l), R1D(l) complex*16 radial function of the 1st kind
c                     R_mn^(1)(c,KSI) and its derivative d/dksi
c       R2(l), R2D(l) complex*16 radial function of the 2nd kind and
c                     its derivative (IKIND = 2), = 0 for IKIND = 1
c       SA(j,l), SAD(j,l)  complex*16 angular function S_mn(c,ETA(j))
c                     and dS_mn/deta (not normalised; divide by
c                     RNORM(l) or sqrt(RNORM(l)) as required)
c       WERR(l)       real*8 |W/W_exact - 1| (Wronskian check) for
c                     IKIND = 2; = -1 for IKIND = 1
c       METH(l)       method used for the radial functions:
c                     11 - expansion in Legendre functions,
c                     13 - Legendre (1st kind) + integration of the
c                          ODE (2nd kind), prolate,
c                     22 - expansion in spherical Bessel functions,
c                     33 - integration of the ODE (oblate),
c                     44 - Jaffe-type expansion (prolate)
c       IER           = 0  o.k.,
c                     < 0  results returned, but one of the internal
c                          routines reported a problem for some
c                          function (check WERR),
c                     > 0  wrong input, nothing computed:
c                          1 KTYPE, 2 M, 3 NE, 4 C, 5 KSI, 6 IKIND,
c                          7 NETA/ETA, 8 LDD/LDS, 9 too many terms
c
c     Notes:
c     - Method selection, accuracy control and all numerical details
c       follow HOMFUNQ of hom6_6q.  One deliberate change: in HOMFUNQ
c       a failed Bessel-type calculation (Wronskian error > EPSW) of a
c       prolate function with n-m+1 <= c falls through to the oblate
c       ODE integrator (label 33); here the Bessel-type result is kept.
c     - Eigenvalues: SEGV1 (Zhang & Jin) for real C, CDLAM_m (Martin,
c       complex tridiagonal eigenproblem) for complex C, as in hom6_6q.
c     - Some internal routines write diagnostics to Fortran unit 7
c       (as in hom6_6q).  Routine names (F02AJF, M01AAF, ...) are those
c       of the old NAG routines bundled with hom6_6q; avoid linking
c       with the NAG library.
c     - Common blocks used: /K1/,/K2/,/FF/,/LEG/,/LEG1/,/FACT/,/PI/,
c       /EPS1/,/EPS3/,/EPS33/,/F2/,/Print/.
c
c=======================================================================

      SUBROUTINE SPHINI
c     initialise the constants in common blocks (cf. main of hom6_6q)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z), COMPLEX*16 (R-S)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /EPS1/ EPS1
      COMMON /EPS3/ EPS3
      COMMON /EPS33/ EPS33
      COMMON /FACT/ FACT(300)
      COMMON /PI/ PI
      COMMON /Print/ iprint
      S=(0D0,0D0)
      S1=(0D0,1D0)
      PI=4d0 * datan(1d0)
      eps3 = 1d-200
      eps33 = eps3
      eps1 = 1d-10
      iprint = 0
c factorial (FACT(I) = (I-1)!, scaled by 1d-300 for I >= 170)
      FACT(1)=1D0
      FACT(2)=1D0
      DO 130 I=3,170
  130 FACT(I)=FACT(I-1)*(I-1D0)
      FACT(170)=FACT(170)*(1.D-300)
      DO 140 J = 170, 297
      FACT(J+1)=J*FACT(J)
  140 CONTINUE
      RETURN
      END

c************************************************
C SPHFUN
C
      SUBROUTINE SPHFUN(KTYPE,M,NE,C,KSI,IKIND,NETA,ETA,EPSW,NKIN,
     *                  RLAM,DCOEF,LDD,RNORM,R1,R1D,R2,R2D,
     *                  SA,SAD,LDS,WERR,METH,IER)
      parameter (nterms=330, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z), COMPLEX*16 (R-S)
      REAL*8 KSI, ETA(*), WERR(*)
      INTEGER METH(*)
      COMPLEX*16 C, RLAM(*), DCOEF(LDD,*), RNORM(*),
     *           R1(*), R1D(*), R2(*), R2D(*),
     *           SA(LDS,*), SAD(LDS,*)
      COMPLEX*16 RDC2(4*nterms), BDC2(4*nterms)
      DIMENSION BDC1(4*nterms)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /FF/ ASINA(nalphas), ACOSA(nalphas)
      COMMON /EPS1/ EPS1
      COMMON /EPS3/ EPS3
      COMMON /EPS33/ EPS33
c
c..... check input
      IER = 0
      IF(KTYPE.NE.0.AND.KTYPE.NE.1) IER = 1
      IF(M.LT.1) IER = 2
      IF(NE.LT.1.OR.NE.GT.nterms) IER = 3
      IF(DREAL(C).LE.0D0) IER = 4
      IF(IKIND.EQ.2.AND.DIMAG(C).NE.0D0) IER = 4
      IF(KTYPE.EQ.0.AND.KSI.LE.1D0) IER = 5
      IF(KTYPE.EQ.1.AND.KSI.LT.0D0) IER = 5
      IF(IKIND.NE.1.AND.IKIND.NE.2) IER = 6
      IF(NETA.LT.0.OR.NETA.GT.nalphas) IER = 7
      IF(IER.NE.0) RETURN
      DO 2 J = 1, NETA
      IF(DABS(ETA(J)).GT.1D0) IER = 7
    2 CONTINUE
      IF(LDS.LT.MAX0(NETA,1)) IER = 8
      IF(IER.NE.0) RETURN
c
c..... constants of /K1/ (as in the main program of hom6_6q)
      K = KTYPE
      IF(K.EQ.0) THEN
        AKSI = KSI**2 - 1D0
        AK = 1D0
      ELSE
        AKSI = KSI**2 + 1D0
        AK = -1D0
      END IF
      EPS1 = EPSW
      IF(EPS1.LE.0D0) EPS1 = 1D-10
      EPS33 = EPS3
c..... number of expansion coefficients (cf. NK in hom6_6q)
      IF(NKIN.GT.0) THEN
        NK = NKIN
      ELSE
        NK = NE + 40
        IF(NK.LT.60) NK = 60
        IF(KSI.GT.1.5D0) NK = NK + 40
      END IF
      IF(NK.GT.nterms-1) NK = nterms - 1
      IF(NK.LE.(NE-1)/2+2) IER = 9
      IF(LDD.LT.NK) IER = 8
      IF(DIMAG(C).NE.0D0.AND.
     *   NE+1+3*IDINT(DSQRT(CDABS(C))).GT.nterms) IER = 9
      IF(IER.NE.0) RETURN
c..... angular arguments
      nal = NETA
      DO 3 J = 1, NETA
      ACOSA(J) = ETA(J)
      ASINA(J) = DSQRT(DABS(1D0 - ETA(J)**2))
    3 CONTINUE
c
c..... eigenvalues lambda_mn(c), n = m, ..., m+ne-1
      IF(DIMAG(C).EQ.0D0) THEN
        KD = 1
        IF(K.EQ.1) KD = -1
        DO 5 JJ = 1, NE
        CALL SEGV1(M,JJ,C,KD,CV,RLAM)
    5   CONTINUE
      ELSE
        IE = 0
        CALL CDLAM_m(RLAM,M,NE,C,IE)
        IF(IE.NE.0) IER = -1
      END IF
c
c..... associated Legendre functions P_{m+i-1}^m(eta) -> /LEG1/
      IF(NETA.GT.0) CALL funlegnn1(M,2*NK+2)
c
c..... choice of the method (as in HOMFUNQ)
      ncc  = 2 * NE / 3
      ncc0 = IDINT(DREAL(C))
      IF(K.EQ.0) THEN
        ifun1 = 44
        IF(KSI.GE.1.3D0) ifun1 = 22
      ELSE
        ifun1 = 33
        IF(KSI.GE.1.5D0) ifun1 = 22
      END IF
      IDC = IKIND
      W1 = 1D0 / (DREAL(C) * (KSI**2 - 1D0 + 2*K))
c     ileg = 1 : Legendre functions of KSI (/LEG/) are computed
      ileg = 0
      DO 6 I = 1, 4*nterms
      BDC1(I) = 0D0
      BDC2(I) = S
    6 CONTINUE
c
c=======================================================================
c     loop over n = m, ..., m+ne-1
c=======================================================================
      DO 1 L = 1, NE
      LM = L + M - 1
      RL = RLAM(L)
      IERR = 0
      IERC = 0
      RA1 = S
      RA2 = S
      RB1 = S
      RB2 = S
      WERR(L) = -1D0
c
c..... expansion coefficients d_r and norm N_mn
      CALL cdcof3(RL,RDC2,NK,BDC1,NK,IDC,RVN,1,M,LM,C,IERC)
      IF(IERC.NE.0) IER = -1
      IF(IDC.EQ.2) THEN
        IF(DABS(BDC1(M+1)).LE.EPS33*1D20) EPS33 = EPS33 / 1D20
      END IF
      DO 7 I = 1, NK
      DCOEF(I,L) = RDC2(I)
    7 CONTINUE
      RNORM(L) = RVN
c
c..... angular functions
      DO 8 J = 1, NETA
      CALL CASWF(SA(J,L),SAD(J,L),M,LM,J,ETA(J),RDC2,NK)
    8 CONTINUE
c
c..... radial functions
      ifl = ifun1
      IF(K.EQ.1.AND.L.GT.ncc.AND.KSI.GT.0.25D0.AND.KSI.LT.0.8D0)
     *   ifl = 11
c
      IF(IKIND.EQ.1) GO TO 600
c
c-----------------------------------------------------------------------
c     IKIND = 2 : 1st and 2nd kind, real c (c1 part of HOMFUNQ)
c-----------------------------------------------------------------------
      W0 = -1D0
      IF(ifl.EQ.11) GO TO 411
      IF(ifl.EQ.22) GO TO 422
      IF(ifl.EQ.33) GO TO 433
      GO TO 444
c
c..... LEG/LEG
  411 CONTINUE
      IF(ileg.EQ.0) CALL funlegnn(M,KSI,NK)
      ileg = 1
      CALL CDRF12(RA1,RA2,AR3,AR4,1,M,LM,C,RDC2,NK,BDC1,NK,IE)
      IF(IE.NE.0) IERR = IE
      RB1 = AR3
      RB2 = AR4
      W = RA1 * AR4 - RA2 * AR3
      W2A = DABS(W / W1 - 1D0)
      WERR(L) = W2A
      METH(L) = 11
      IF(W0.GE.0D0.AND.W2A.GT.W0.AND.K.EQ.0.AND.ifl.EQ.44) GO TO 413
      GO TO 500
c
c..... LEG/DIFF (prolate): 2nd kind by integration of the ODE
  413 CONTINUE
      CALL cDRSF20(RR3,RR4,KSI,RA1,RA2,M,C,RL,IE)
      IF(IE.NE.0) IERR = IE
      W = RA1 * RR4 - RA2 * RR3
      W2B = DABS(W / W1 - 1D0)
      IF(W2B.LT.W2A) THEN
        RB1 = RR3
        RB2 = RR4
        WERR(L) = W2B
        METH(L) = 13
      END IF
      GO TO 500
c
c..... BES/BES
  422 CONTINUE
      CALL CDRB12aa(RA1,RA2,AR3,AR4,1,M,LM,C,KSI,RDC2,NK,IE)
      IF(IE.NE.0) IERR = IE
      RB1 = AR3
      RB2 = AR4
      W = RA1 * AR4 - RA2 * AR3
      W2A = DABS(W / W1 - 1D0)
      WERR(L) = W2A
      METH(L) = 22
      IF(W2A.LE.EPS1) GO TO 500
      IF(K.EQ.0) THEN
        IF(L.GT.ncc0) GO TO 411
        GO TO 500
      END IF
c
c..... DIFF/DIFF (oblate): integration of the ODE
  433 CONTINUE
      DO 9 I = 1, M+1
      BDC2(I) = BDC1(I)
    9 CONTINUE
      CALL cDRSF212(RA1,RA2,RR3,RR4,2,M,LM,C,KSI,RL,
     *              RDC2,NK,BDC2,M+1,IE)
      IF(IE.NE.0) IERR = IE
      RB1 = RR3
      RB2 = RR4
      W = RA1 * RR4 - RA2 * RR3
      W2A = DABS(W / W1 - 1D0)
      WERR(L) = W2A
      METH(L) = 33
      IF(W2A.GT.EPS1.AND.L.GT.ncc0) GO TO 411
      GO TO 500
c
c..... JAF/JAF (prolate)
  444 CONTINUE
      CALL CDRG1bn(RA1,RA2,AR3,AR4,1,M,LM,C,KSI,
     *             RDC2,NK,BDC1,NK,RL,IE)
      IF(IE.NE.0) IERR = IE
      RB1 = AR3
      RB2 = AR4
      W = RA1 * AR4 - RA2 * AR3
      W0 = DABS(W / W1 - 1D0)
      WERR(L) = W0
      METH(L) = 44
      IF(W0.GT.EPS1.AND.L.GT.ncc0) GO TO 411
      GO TO 500
c
c-----------------------------------------------------------------------
c     IKIND = 1 : 1st kind only, real or complex c (c2 part of HOMFUNQ)
c-----------------------------------------------------------------------
  600 CONTINUE
      ifun2 = 1
      IF(ifl.EQ.22.OR.KSI.GT.1.6D0) ifun2 = 2
      IF(ifl.EQ.33) ifun2 = 2
      IF(ifl.EQ.44) ifun2 = 3
      GO TO (731, 732, 733), ifun2
c..... Legendre functions
  731 CONTINUE
      IF(ileg.EQ.0) CALL funlegnn(M,KSI,NK)
      ileg = 1
      CALL CDRF12(RA1,RA2,V,V,0,M,LM,C,RDC2,NK,BDC1,NK,IE)
      METH(L) = 11
      GO TO 740
  732 CONTINUE
      IF(ifl.EQ.22.OR.KSI.GT.1.6D0) THEN
c..... spherical Bessel functions
        CALL CDRB12aa(RA1,RA2,V,V,0,M,LM,C,KSI,RDC2,NK,IE)
        METH(L) = 22
      ELSE
c..... integration of the ODE (oblate)
        CALL cDRSF212(RA1,RA2,RRV,RRV,0,M,LM,C,KSI,RL,
     *                RDC2,NK,BDC2,M+1,IE)
        METH(L) = 33
      END IF
      GO TO 740
c..... Jaffe-type expansion (prolate)
  733 CONTINUE
      CALL CDRG1bn(RA1,RA2,AR3,AR4,0,M,LM,C,KSI,
     *             RDC2,NK,BDC1,NK,RL,IE)
      METH(L) = 44
  740 CONTINUE
      IF(IE.NE.0) IERR = IE
c
  500 CONTINUE
      R1(L) = RA1
      R1D(L) = RA2
      R2(L) = RB1
      R2D(L) = RB2
      IF(IERR.NE.0) IER = -1
    1 CONTINUE
      RETURN
      END

c************************************************
C CASWF
C   angular spheroidal function S_mn(c,eta) and dS_mn/deta for
C   complex expansion coefficients AD (cf. DASWF1 of hom6_6q);
C   P_{m+i-1}^m(eta) are taken from /LEG1/ (computed by funlegnn1)
C
      SUBROUTINE CASWF(RS,RD,M,N,ICOS,X,AD,IA)
      parameter (nterms=330, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z), COMPLEX*16 (R-S)
      COMPLEX*16 AD(IA)
      COMMON /LEG1/ P(nalphas,2*nterms), PD(nalphas,2*nterms)
      eps4 = 1d-40
      RS = (0D0,0D0)
      RD = (0D0,0D0)
      IF(1D0-DABS(X).LT.1D-12) RETURN
      JK = MOD(N-M,2)
      DO 13 I = 1, IA
      J1 = 2*(I-1) + 1 + JK
      RSI = AD(I) * P(ICOS,J1)
      RDI = AD(I) * PD(ICOS,J1)
      RS = RS + RSI
      RD = RD + RDI
      IF(I.GT.2.AND.CDABS(RSI).LE.eps4*CDABS(RS).AND.
     *   CDABS(RDI).LE.eps4*CDABS(RD)) GO TO 14
   13 CONTINUE
   14 CONTINUE
      RETURN
      END

c=======================================================================
c     Routines below are copied without changes from hom6_6q
c=======================================================================

c************************************************
C funlegnn1
C
      SUBROUTINE funlegnn1(M,NF)
      parameter (nterms=330,nalphas=19)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 s, s1
      DIMENSION P(nalphas,2*nterms),PD(nalphas,2*nterms),
     *          Pp(2*nterms),PDp(2*nterms)
      COMMON /LEG1/ P,PD
      COMMON /K1/ S, S1, AKSI, AK, kkK, NK, nal
      COMMON /FF/ ASIN(nalphas),ACOS(nalphas)
      IF(NF.GT.2*nterms) WRITE(7,1) NF
      IF(NF.GT.2*nterms) WRITE(*,1) NF
    1 FORMAT(1X,'funlegnn1 NF > 250, NF=',2I10)
c      IF(NF.GT.2*nterms) GO TO 900
      IF(NF.GT.2*nterms) stop 666
c      IF(NF.GT.250) nf = 250
         DO j  = 1, nalphas
         DO I  = 1, 2*nterms
         P (j,I) = 0D0
         PD(j,I) = 0D0
         end do
         end do
      do 6 ia=1,nal
      x=acos(ia)
c      IF(M.EQ.1) THEN
c        CALL DUGOL1(X,Nf,Pp,PDp)
c      ELSE
        CALL DUGOLM1(X,M,Nf,Pp,PDp)
c      ENDIF
        A1=DSQRT(1.D0-X*X)
      DO I = 1, Nf
        P(ia,I)=Pp(I)*A1**M
        PD(ia,I)=-PDp(I)*A1**(M-2)
      ENDDO
c      CALL DLEGF1(P(ia,1),D,X,D,d,M,M,0,0,IER1)
c       IF((IER1).NE.0) WRITE(7,2) M,IER1
c    2 FORMAT(1X,'funlegnn1 M,I,IER1,2,3,4=',6I10)
c       DO 3 I=2,NF
c       if(i.gt.2) d=p(ia,i-2)
c      CALL DLEGF1(P(ia,I),PD(ia,I-1),X,P(ia,I-1),d,M+I-1,M,0,1,IER3)
c
c       print *,nf,i,P(ia,I)
c
c       IF((IER3).NE.0) WRITE(7,2) M,I,IER3
c    3 CONTINUE
    6 CONTINUE
  900 RETURN
      END

c************************************************
C cdcof3
C
      SUBROUTINE cdcof3(VL,AD,IA,BD,IB,ID,VN,IN,M,N,C,IER)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 VL,AD,VN,C
      DIMENSION AD(IA),BD(IB),RL(4*nterms)
      COMMON /EPS3/ EPS3
      COMMON /EPS33/ EPS33
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      common /FACT/ FACT(300)
      AC=C
      IER=0
c      eps33 = eps3
c+&      eps4 = 1d-20
               eps4 = 1d-40
C* ЏPOBEPKA BBO„€M›X ‡HA—EH€‰ ЏAPAMETPOB
      IF(M.LT.0.OR.N.LT.M) IER=1
      IF(AC.LT.0D0) IER=3
      IF(ID.EQ.1.AND.IA.LE.((N-M)/2+2)) IER = 4
      IF(ID.EQ.0.AND.IN.NE.0) IER = 6
      IF(IER.NE.0) GO TO 900
      aM = M
      aMM = 2d0 * M
      S2=(1D0,0D0)
      bM = 4D0 * am**2 - 1D0
      RC2=C*C
      RC4=RC2*RC2
      N1=N-M
      IF(ID.EQ.0) GO TO 900
C**** B›—€C‹EH€E K-TOB PA‡‹O†EH€џ
      IF(MOD(N1+1,2)) 21,20,21
   20 I11=1
      GO TO 22
   21 I11=0
   22 CONTINUE
      L1=(N1-I11)/2
      L11=L1+1
      RL(L11)=S2
      IF(L1.EQ.0) GO TO 24
C***  B›—€C‹EH€E K-TOB ЏP€ R<N-M
      DO 23 J=1,L1
      JJ1=N1-2*J
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(amm+2D0*JJ1-1D0)/
     *   (amm+2D0*JJ1+3D0))        -VL
      RV=S2/RB
      RW=RV
      JJ2=JJ1
      IF(JJ2.LE.1) GO TO 103
  102 CONTINUE
      II=2*M+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2-2D0)*(am+JJ2-1D0)+AK*RC2/2D0*(1D0-bm/(II1-5D0)/
     *   (II1-1D0))        -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
      IF(CDABS(RV/rw).LE.eps4.and.j.ge.l1/2) GO TO 103
c      IF(CDABS(RV/rw).LE.eps4) GO TO 103
      JJ2=JJ2-2
      IF(JJ2.LE.1) GO TO 103
      GO TO 102
  103 RL(L11-J)=-(amm+JJ1+1D0)*(amm+JJ1+2D0)*AK*RC2/(amm+2D0*JJ1+3D0)/
     *           (amm+2D0*JJ1+5D0)*RW
   23 RL(L11-J)=RL(L11-J)*RL(L11+1-J)
   24 CONTINUE
C***  B›—€C‹EH€E K-TOB ЏP€ R>N-M
      IA1=IA-1
      DO 25 J=L11,IA1
      JJ1=N1+2*(J-L11)+2
      II=2*M+JJ1
      II1=II+JJ1
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(II1-1D0)/(II1+3D0))
     *   -VL
      RV=-JJ1*(JJ1-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RW=RV
      JJ2=JJ1
  104 CONTINUE
      JJ2=JJ2+2
      II=2*m+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2)*(am+JJ2+1D0)+AK*RC2/2D0*(1D0-bm/(II1-1D0)/(II1+3D0))
     *          -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV

c        if(n.eq.1.and.dimag(c).gt.0.1d0.and.j.le.3)  write(7,*) j, rw

c       IF(CDABS(RV/rw).LE.eps4.and.j.ge.ia1/2) GO TO 105
        IF(CDABS(RV/rw).LE.eps4) GO TO 105
c       IF(CDABS(RV/rw).LE.eps3) GO TO 105
      GO TO 104
  105 RL(J+1)=RW*(amm+2D0*JJ1-1D0)*(amm+2D0*JJ1+1D0)/((amm+JJ1-1D0)*
     *      (amm+JJ1)*AK*RC2)
      RL(J+1)=RL(J)*RL(J+1)
       IF(J.GT.L11) RL(J+1)=RL(J+1)/RL(L11)
       IF(CDABS(RL(J+1)).LT.eps3) GO TO  125
   25 CONTINUE
      GO TO 126
  125 IJ=J+2
         DO 127 I=IJ,IA
  127    RL(I)=S
  126 CONTINUE

C *                  B›—€C‹EH€E K-TA D(N-M)
      L4=(N+M+I11)/2
c      if(l1+1.gt.170) stop 1701
      F2=1D0/Fact(l1+1)
      DO  I = 1, L4
       F2 = F2 * (L4 + I)
       if(dabs(f2).gt.1d300) go to 27
      end do
      iff = 0
      go to 270
   27 continue
      f2 = f2 * 1d-300
      iff = 1
c Draine!
c      DO Ij = i, l4
      DO Ij = i+1, l4
c       F2 = F2 * (J3 + Ij)
       F2 = F2 * (L4 + Ij)
      end do
 270  continue
      il1 = 1
      if(mod(l1,2).ne.0) il1 = -1
      D1 = iL1 * F2 / 2D0**N1
       if(iff.eq.1) d1 = d1 * 1d300
      SNM=S

      DO 28 J=1,IA
      J1=2*(J-1)+I11
      J2=(J1-I11)/2
      J3=(J1+2*M+I11)/2
      rF2=RL(J)/2D0**J1
      DO I=1,J3
       rF2 = rF2 * (J3 + I)
       if(cdabs(rf2).gt.1d300) go to 30
      end do
      iff = 0
      go to 330
   30 continue
      rf2 = rf2 * 1d-300
      iff = 1
c DRAINE+!
c      DO Ij = i, J3
      DO Ij = i+1, J3
       rF2 = rF2 * (J3 + Ij)
      end do
 330  continue
      rF1=rF2
      DO 29 I=1,J2
   29 rF1=rF1/I
       SG=(-1D0)**J2*rF1
       if(iff.eq.1) SG = sg * 1d300
       IF(CDABS(SG).LT.eps3.AND.J.GT.ia/2) GO TO 128
   28 SNM=SNM+SG
  128 AD(L11)=D1/SNM

      DO 31 I=1,IA
   31 AD(I)=RL(I)*AD(L11)
      IF(ID.NE.2) GO TO 901
C*       B›—€C‹EH€E KOќ””€–€EHTOB PA‡‹O†EH€џ ЏP€  R<0
      BD(1)=AD(1)
      DO 223 J=2,IB
      JJ1=I11-2*(J-1)
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(amm+2D0*JJ1-1D0)/
     *   (amm+2D0*JJ1+3D0))        -VL
      RV=1D0/RB
      RW=RV
      JJ2=JJ1
  224 CONTINUE
      II=2*M+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2-2D0)*(am+JJ2-1D0)+AK*RC2/2D0*(1D0-bm/(II1-5D0)/
     *   (II1-1D0))        -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
c          print *, j, RV/rw
c          pause
c      IF(CDABS(RV/rw).LE.eps4.and.j.ge.ib/2) GO TO 225
      IF(CDABS(RV/rw).LE.eps4) GO TO 225
      JJ2=JJ2-2
      GO TO 224
  225 CONTINUE
      IF(JJ1.EQ.-2*M-2.OR.JJ1.EQ.-2*M-1) GO TO 226
      BD(J)=-(amm+JJ1+1D0)*(amm+JJ1+2D0)*AK*RC2/(amm+2D0*JJ1+3D0)/
     *           (amm+2D0*JJ1+5D0)*RW
      GO TO 227
  226 BD(J)= AK*RC2/(amm+2D0*JJ1+3D0)/(amm+2D0*JJ1+5D0)*RW
      IF(I11.EQ.1) BD(J)=-BD(J)
  227 BD(J)=BD(J)*BD(J-1)
          IF(DABS(BD(J)).LT.eps33) GO TO 1223
  223 CONTINUE
        GO TO 901
 1223 IJ=J+1
        DO 1224 I=IJ,IB
 1224   BD(I)=0D0
  901 IF(IN.EQ.0) GO TO 900
C*****  B›—€C‹EH€E HOPM€PYћ™EѓO MHO†€TE‹џ *
      VN=S
         iff = 0
         if(cdabs(ad(1)).gt.1d150) iff = 1

      DO 32 J=1,IA
      JJ1=2*(J-1)+I11
         rf1 = 1d0
         if(iff.eq.1) rf1 = 1d-300
      rF1 = rf1 * AD(J) / (amm+2D0*JJ1+1D0) * ad(j)
      IF(M.EQ.0) GO TO 333
      DO 33 I=1,2*M
   33 rF1=rF1*(JJ1+I)
  333 RVN=rF1
      VN=VN+RVN
      IF(CDABS(RVN/(vn+1d-50)).LT.eps4.and.j.gt.ia/2) GO TO 34
c      IF(CDABS(RVN/(vn+1d-50)).LT.eps4) GO TO 34
   32 continue
   34 VN=2D0*VN
         if(iff.eq.1) vn = vn * 1d300
  900 RETURN
      END

c************************************************
c cdrb12aa
c
      SUBROUTINE CDRB12aa(R,DR,AR,ADR,ID,M,N,C,X,AD,IA,IER)
c      parameter (nterm=400)
      parameter (nterm=500)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AD,C,S1,BESJ,BESY,AC,SIG,RK,R,DR,F2,RS,RD,S,BR,ABR
      COMPLEX*16 sS1
      DIMENSION AD(IA),BESJ(4*nterm),BESY(4*nterm)
      DIMENSION iiy(4*nterm)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /pi/ pi
      factor=1d250
      IER=0
      IF(X.LT.1D-10) IER=10
      IF(IER.NE.0) GO TO 900
      M2=2*M
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 5
    2 JK=0
    5 CONTINUE
      AM=M
      AM=AM/2D0
      PI2 = pi / 2d0
      X2=X*x
      AC=C*X
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
      F2=AD(I)
      DO 8 I1=1,M2
    8 F2=F2*(J+I1)
    7 SIG=SIG+F2
      RK=CDSQRT(PI2/AC)/SIG*((X2-AK)/X2)**AM
       NUM=3*CDABS(AC)+10
c       if(x.gt.2d0) NUM = 3 * CDABS(AC) + 10 * (m + x)
       if(x.gt.2d0) NUM = 3 * CDABS(AC) + 6 * (m + x)
      if(NUM.lt.60) num=60
      if(NUM-1.ge.4*nterm) write(*,*) 'cdrb12aa ',num, 4*nterm
      if(NUM-1.ge.4*nterm) num = 4*nterm
      IF(ID.EQ.1) CALL CESSEL0(1D0/AC,NUM,BESJ,BESY,iiy)
      IF(ID.EQ.0) CALL BESSJJ0(1D0/AC,NUM,BESJ)
      R=S
      DR=S
      AR=0D0
      ADR=0D0
      F1=AK/X*M/(X2-AK)
      NM2=MIN0(IA,NUM/2-2)

      DO 13 I=1,NM2
      J=2*(I-1)+JK
c      F2=AD(I)*S1**(J+M-N)
      F2=AD(I)*sS1(J+M-N)
      DO 18 I1=1,M2
   18 F2=F2*(J+I1)
      IF(I.EQ.1) GO TO 21
c+&      IF(CDABS(RD).LT.eps3.AND.CDABS(RS).LT.eps3.AND.I.GT.(M+5)
c+&     *   .AND.ID.EQ.0)   GO TO 19
   21 RS=F2*BESJ(J+M+1)
      RD=F2*((F1+(M+J)/X)*BESJ(J+M+1)-C*BESJ(J+M+2))
      R=R+RS
      DR=DR+RD
   20 CONTINUE
      IF(ID.EQ.0) GO TO 13
      IF(I.EQ.1) GO TO 22
c+&      IF(CDABS(BR).LT.eps3.AND.CDABS(ABR).LT.eps3.AND.I.GT.(M+5))
c+&     *   GO TO 14
   22 BR = BESY(J+M+1) * F2 * factor**(iiy(j+m+1))
                     AR  = AR  + BR
      ABR=(F1 +(M+J)/X)*BESY(J+M+1) * f2* factor**(iiy(j+m+1))
     *    -C*BESY(J+M+2) * f2 * factor**(iiy(j+m+1))
                     ADR = ADR + ABR
   13 CONTINUE
   14 AR= AR* RK
      ADR= ADR* RK
   19 R=R*RK
      DR=DR*RK

c               IF(Id.EQ.0) then
c                       write(*,*) r,rk
c                 pause
c               end if
  900 RETURN
      END

c************************************************
c CDRF12
c
      SUBROUTINE CDRF12(R,DR,AR,ADR,I12,M,N,C,AD,IA,BD,IB,IER)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 AD,C,DR
      DIMENSION AD(IA ),BD(IB+1)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /LEG/ P(2*nterms),PD(2*nterms),Q(2*nterms),QD(2*nterms)
      COMMON /FACT/ FACT(300)
      COMMON /eps3/ eps3
       eps4 = 1d-20
      ier=0
      M2=2*M
      aM2 = 2d0 * M
      nm1=n-m
      nm2=n+m
      IF(MOD((N-M),2)) 1,2,1
    1  CONTINUE
      JK=1
      GO TO 36
    2  CONTINUE
      JK=0
   36 CONTINUE
      IF(MOD(N,2)) 3,4,3
    3  CONTINUE
      Jn=1
      GO TO 5
    4  CONTINUE
      Jn=0
    5 CONTINUE
      IF(MOD(m,2)) 33,34,33
   33  CONTINUE
      Jm=1
      GO TO 35
   34  CONTINUE
      Jm=0
   35 CONTINUE
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
       rf2=AD(I)*fact2(m2+j,j)
      SIG=SIG+rf2
      if(cdabs(rf2/sig).lt.eps4) go to 8
    7  CONTINUE
    8  CONTINUE
      if((nm1-jk)/2+1.gt.170) stop 1702
      ff3 =  fact3(nm2+jk,(nm2+jk)/2,1d0/fact((nm1-jk)/2+1))
      RK1=(am2+1d0+jk*2d0)/2d0**nm2 * ff3 /
     *    AD(1)/fact(m+1)*sig/C**(M+JK)
      if(k.eq.1) rk1=rk1/(-s1)**(M+JK)*(-s1)**jn
      R=S
      DR=S
      IF(I12.EQ.0) GO TO 70
      RK2=2d0**nm1/(am2-1d0)*fact(m2+1)/
     *    fact(m+1) / ff3
     *    *bd(m+1)*sig/C**(M-1-JK)
      if(jk.eq.1.and.k.eq.0) rk2=-rk2/(am2-3d0)
      if(k.eq.0) go to 339
      if(jk.eq.1.and.jm.eq.0) rk2=rk2/(am2-3d0)
      if(jk.eq.1.and.jm.eq.1) rk2=-rk2/(am2-3d0)*(-s1)
      if(jk.eq.0.and.jm.eq.0) rk2=rk2*s1
      if(k.eq.1) rk2=rk2/(-s1)**(M-1-JK)
  339 continue
      AR=0D0
      ADR=0D0
      DO 12 I=1,M
        J2=2*(M-I)+1+JK
      ADR=ADR+BD(I+1)*QD(J2)
   12 AR=AR+BD(I+1)*Q(J2)           
   70 CONTINUE
      DO 13 I=1,IA
        J1=2*(I-1)+1+JK
        J2=2*(M+I)-1+JK
      RS=AD(I)*P(J1)
      RD=AD(I)*PD(J1)
      R=R+RS
      DR=DR+RD
c      IF(I12.EQ.0) write(7,*) i,r,rs
c      IF(I12.EQ.0) write(7,*) j1,AD(I),P(J1)
      IF(I12.EQ.0) GO TO 17
      BR=AD(I)*Q(J2)

c         write(7,*) m,n,i,AD(I),j2,Q(J2)

      ABR=AD(I)*QD(J2)
      AR=AR+BR
      ADR=ADR+ABR                                  
c????      
c      IF(DABS(BR/ar).LT.eps4.AND.DABS(ABR/adr).LT.eps3) GO TO 14
c+&       IF(DABS(BR/ar).LT.eps4.AND.DABS(ABR/adr).LT.eps4) GO TO 14
   17 CONTINUE
      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RD/dr).LT.eps4) GO TO 14
c+&      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RD/dr).LT.eps4) GO TO 14
   13 CONTINUE
   14 CONTINUE
      IF(I12.EQ.0) GO TO 71
      MR=M+1
      DO 15 I=MR,IB
        J2=2*(I-M)-JK
      if (k.eq.1.and.mod(m,2).eq.0) go to 116
      BR=BD(I+1)*P(J2)

c         write(7,*) m,n,i+1,bd(I+1),j2,P(J2)

      ABR=BD(I+1)*PD(J2)
      go to 115
  116 CONTINUE
      BR=-BD(I+1)*P(J2)

c         write(7,*) m,n,i+1,bd(I+1),j2,P(J2)

      ABR=-BD(I+1)*PD(J2)
  115 CONTINUE
      AR=AR+BR
c      if(m.eq.29) then
c      write(*,*) m,n,i+1,bd(I+1),j2,P(J2)
c      write(*,*) i,ar,br      
c      end if
      ADR=ADR+ABR
c+&      IF(DABS(BR/ar).LT.eps4.AND.DABS(ABR/adr).LT.eps4) GO TO 16
   15 CONTINUE
   16 AR=AR/RK2
c        if(m.eq.29) then
c          write(*,*) i,ar,br,rk2
c        end if
      ADR=ADR/RK2
   71 R=R/RK1

c               IF(I12.EQ.0) then
c                       write(*,*) r,rk1
c                 pause
c               end if

      DR=DR/RK1
  900 RETURN
      END

c************************************************
C DLEGf1
C   
      subroutine DLEGF1(rs,rd,X1,x2,x3,n,M,ic,id,iER)
      IMPLICIT REAL*8(A-H,O-Z) 
      complex*16 s, s1
      COMMON /FACT/ FACT(300)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
      ier=0
      res=0d0
      red=0d0
      eps4 = 1d-20
      xx2=dabs(x2)
      nm=n-m
      j=mod(n,4)
      jm=mod(m,4)
      if(n.eq.0.and.m.eq.0) res=1d0
      if(n.eq.0.and.m.eq.0)  go to 100
      if(n.le.0.or.m.lt.0) ier=1
      if(id.ne.0.and.dabs(x1-1d0).lt.1d-6.and.ic.eq.0) ier=2
      if(ier.ne.0) go to 100
      if(dabs(x1).gt.1d-6)  go to 20
c**************************************   << x1 = 0 >>
      IF(MOD(iabs(nm),2)) 100,1,100
    1 continue
      k1=(m+n)/2
      e=1d0
      do 2 k=1,k1
    2 e=e*(k1+k)
      if(nm/2+1.gt.170) stop 1703
      res=e/fact(nm/2+1)/2d0**n
      if(id) 50,102,50
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   << x1 /= 0 >>
   20 continue
      if(m.gt.n) go to 100
      if(nm.ge.2.and.x1.gt.1d0.and.ic.eq.0) go to 48
      z=1d0/(2d0**n)
      k1=(nm+2)/2
      d=0d0
      if(ic.ne.0) go to 21
c.....................................   << ic = 0 >>
c      if(nm.gt.169) go to 291
      if(nm.gt.169) stop 170

c       print *, m, n, k1
c       pause

      do 22 k=1,k1
      k3=nm-2*(k-1)
      if(k3+1.gt.170) stop 1704
      h=1d0/fact(k3+1)
      k4=n-k+1

c       print *,  k, k4
c       pause

      do 27 i=1,k4
      if(dabs(h).gt.1d300) print *,  n, h

   27 h=h*(k4+i)
      id = 1
       if(mod(k-1,2).ne.0) id = -1
       dd = id * h / fact(k) * x1**k3
c       print *,  dd, d
      d = d + dd
c+&      if(dabs(dd/d).lt.eps4) go to 23
   22 continue
   23 continue
       res=d*z*dsqrt(dabs(1d0-x1**2)**m)
      if(id) 50,102,50     
c<<<<<<<<<++++++++++++++++++++++++++++   << ic /= 0 >>
   21 continue
c      if(nm.gt.169) go to 29
      if(nm.gt.169) stop 1705
      do 28 k=1,k1
      k3=nm-2*(k-1)
      h=1d0/fact(k3+1)
      k4=n-k+1
c      if(k4.gt.169) stop 17055
      do 33 i=1,k4
   33 h=h*(k4+i)
   28 d=d+h/fact(k)*x1**k3
      res=d*z*dsqrt((1+x1**2)**m)
      if(id) 50,102,50     
   48 continue
      res=((2d0*n-1d0)*x1*x2-(n+m-1d0)*x3)/nm
      if(id) 50,102,50     
   50 continue
c=========================             <<  derivature  >>
c*****************************************   << x1 /= + - 1 >>
      if(ic.ne.0) go to 101
      red=((nm)*res-x1*n*x2)/(x1**2-1d0)
       rs=res
       rd=red
       return
  101 continue
      red=((nm)*res-x1*n*xx2)/(x1**2-ak)
  100 CONTINUE        
      if(jm.eq.2.or.jm.eq.3) go to 103
      rs=res
      if(j.eq.2.or.j.eq.3) rs=-dabs(res)
      go to 104
  103 continue
      rs=-dabs(res)
      if(j.eq.0.or.j.eq.1) rs=-rs
  104 continue
      j=mod(n-1,4)                          
      rd=red
      if(j.eq.2.or.j.eq.3) rd=-dabs(red)
      RETURN      
  102 continue
       rs=res
      if(ic.eq.0.and.jm.eq.3.and.x1.gt.1d0) rs=dabs(res)
      if(ic.ne.0.and.jm.eq.2) rs=-dabs(res)
      if(ic.ne.0.and.jm.eq.3) rs=-dabs(res)
      RETURN                                               
      END

c************************************************
C DLEGf2
C   
      subroutine DLEGF2(res,red,X1,x2,n,M,ic,id,iER) 
      IMPLICIT REAL*8(A-H,O-Z)  
      complex*16 s1, s
c + &
      dimension r(3000),q(3000)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
      COMMON /FACT/ FACT(300)
      COMMON /PI/ PI
      ier=0
      res=0d0
      red=0d0
      if(n+m.lt.0) ier=1
      if(ic.eq.0.and.dabs(x1).le.1d0) ier=2
      if(ic.ne.0.and.dabs(x1).lt.0d0) ier=2
      if(n+m-1.lt.0.and.id.ne.0) ier=3
      if(ier.ne.0) go to 100
      if(dabs(x1).gt.1d-6.or.ic.eq.0)  go to 20
c**************************************   << x1 = 0 & ic /= 0 >>
      if(m.gt.n) go to 3
c                             m  < = n
      k1=n-m
      k2=n+m-1
      d=1d0
      e=1d0
      IF(MOD(k1,2)) 2,1,2
    1 continue
      if(k1.eq.0) go to 6
      do 4 i=1,k2,2
    4 d=d*i
      do 5 i=2,k1,2
    5 e=e*i
    6 res=-pi/2d0*d/e
      if(id) 50,100,50
    2 continue
      do 7 i=2,k2,2
    7 d=d*i
      do 8 i=1,k1,2
    8 e=e*i
      res=d/e
      if(id) 50,100,50
c                                m  >  n 
    3 continue
      IF(MOD(iabs(N-M),2)) 10,9,10
    9 go to 100
   10 continue
      k1=(n-m-1)/2
      k2=(m+n-1)/2
      d=1d0
      do 11 i=1,k1
   11 d=d*(k1+i)
      if(k1+1.gt.170) stop 1706
      res=2d0**n*d*fact(k1+1)
      if(id) 50,100,50
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   << x1 /= 0 >>
   20 continue
      if(m.eq.0)  go to 61
      if(m.eq.1)  go to 71
c *********************************     m  > 1
      if(ic.eq.0) res=dlog((x1+1d0)/(x1-1d0))/2d0
      if(ic.ne.0.and.x1.lt.0.5d0) res=(datan(x1)-pi/2d0)
      if(ic.ne.0.and.x1.ge.0.5d0) res=-datan(1d0/x1)
      r(1)=-ak/(dsqrt(x1**2-ak))
      r(2)=-2d0*x1*ak/(dsqrt(x1**2-ak))*r(1)
      if(m.eq.2) go to 221
      do 21 i=3,m
   21 r(i)=-2d0*(i-1d0)*x1*ak/(dsqrt(x1**2-ak))*r(i-1)+
     *     (2d0-i)*(i-1d0)*r(i-2)
  221 continue
      q0=r(m)
      if(n.eq.0) res=q0
      if(n.eq.0) go to 26
      j=n
      if(n.lt.0) j=-n
      nmax=30*x1+j+50
c+&      if(nmax.gt.500) nmax=500
c      if(nmax.gt.1500) print *, 'dlegf2', nmax
c      if(nmax.gt.1500) nmax=1500
      if(nmax.gt.3000) print *, 'dlegf2', nmax
      if(nmax.gt.3000) nmax=3000
      r(nmax)=x1-dsqrt(x1**2-ak)
      xer=dabs(x1-1d0/dsqrt(3d0))
      do 22 i=nmax-1,1,-1
   22 r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
      q(1)=r(1)*q0*ak
      if(xer.lt.1d-10.and.m.eq.3) q(1)=3d0
      if(xer.lt.1d-10.and.m.eq.6) q(1)=-144d0
      if(xer.lt.1d-10.and.m.eq.9) q(1) = 4.536d4
      if(xer.lt.1d-10.and.m.eq.12) q(1) = -4.35456d+07
      if(xer.lt.1d-10.and.m.eq.15) q(1) = 9.340529643126767d+10
      if(xer.lt.1d-10.and.m.eq.18) q(1) = -3.76610155207262d+14
      if(xer.lt.1d-10.and.m.eq.21) q(1) = 2.554546682742534d+18
      if(xer.lt.1d-10.and.m.eq.24) q(1) = -2.697601737669382d+022
      if(xer.lt.1d-10.and.m.eq.27) q(1) =  4.188026680200701d+026
      if(xer.lt.1d-10.and.m.eq.30) q(1) = -9.146650294330210d+030
      if(xer.lt.1d-10.and.m.eq.33) q(1) =  2.713536749598336d+035
      if(xer.lt.1d-10.and.m.eq.36) q(1) = -1.062838075962123d+040
      if(xer.lt.1d-10.and.m.eq.39) q(1) =  5.367863692727776d+044
      if(xer.lt.1d-10.and.m.eq.42) stop 42
      if(xer.lt.1d-10.and.mod(m,3).eq.0) q0=0d0
      if(n.eq.1) go to 24
      if(n.lt.0) go to 225
      do 23 i=1,j-1
   23 q(i+1)=r(i+1)*q(i)*ak**(i+1)
       go to 24
c******************************************      n < 0
  225 continue 
      if(n.eq.-1) res=(ak*x1*q0+(m-1d0)*q(1))/m
      if(n.eq.-1) go to 26
      r(1)=q0
      r(2)=(ak*x1*q0+(m-1d0)*q(1))/m

c         print *, m,n,q0,q(1)
c         write(7,*) m,q(1)

      do 25 i=-2,n,-1
      r(1-i)=(ak**i*(2d0*i+3d0)*x1*r(-i)-
     *       (i-m+2d0)*r(-i-1))/(i+m+1d0)
   25 continue 
      res=r(-n+1)

c         print *,res

      go to 26
   24 continue
      res=q(n)
   26 continue
      if(id) 50,100,50
c***********************************    m = 0
   61 continue
      if(ic.eq.0) res=dlog((x1+1d0)/(x1-1d0))/2d0
      if(ic.ne.0.and.x1.lt.0.5d0) res=(datan(x1)-pi/2d0)
      if(ic.ne.0.and.x1.ge.0.5d0) res=-datan(1d0/x1)
      if(n.eq.0) go to 65
      q0=res
      nmax=1.5*x1+n+10
      r(nmax)=x1-dsqrt(x1**2-ak)
      do 62 i=nmax-1,1,-1
      r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
   62 continue
      q(1)=r(1)*q0*ak
      if(n.eq.1) GO TO 64
      do 63 i=1,n-1
      q(i+1)=r(i+1)*q(i)*ak**(i+1)
   63 continue
   64 continue
      res=q(n)
   65 if(id) 50,100,50
c***********************************    m = 1
   71 continue
      res=-ak/(dsqrt(x1**2-ak))
      if(n.eq.0) go to 75
      q0=res
      if(n.eq.-1) res=ak*x1*q0
      if(n.eq.-1) go to 75
      nmax=30*x1+n+50
c+&
      if(nmax.gt.3000) print *, 'dlegf2', nmax
      if(nmax.gt.3000) nmax=3000
      r(nmax)=x1-dsqrt(x1**2-ak)
      do 72 i=nmax-1,1,-1
      r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
   72 continue
      q(1)=r(1)*q0*ak
      if(n.eq.1) GO TO 74
      do 73 i=1,n-1
   73 q(i+1)=r(i+1)*q(i)*ak**(i+1)
   74 continue
      res=q(n)
   75 if(id) 50,100,50
   50 continue
c=========================             <<     Їа®Ё§ў®¤­ п    >>
      if(id.eq.0) go to 100
      if(ic.eq.0.or.mod(n,2).eq.0)  go to 174
      red=-((n-m)*res+x1*n*x2)/(x1**2-ak)
      go to 100
  174 continue
      red=((n-m)*res-x1*n*x2)/(x1**2-ak)
  100 CONTINUE                                                      
      RETURN                                                        
      END

c************************************************
C funlegnn
C
      SUBROUTINE funlegnn(M,X,NF)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 s, s1
      DIMENSION P(2*nterms),PD(2*nterms),Q(2*nterms),QD(2*nterms)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /LEG/ P, PD, Q, QD
      IF(NF.GT.2*nterms) WRITE(7,1) NF
    1 FORMAT(1X,'funlegnn  NF > 2*nterms, NF=',2I10)
      IF(NF.GT.2*nterms) GO TO 900
         DO I  = 1, 2*nterms
         P (I) = 0D0
         PD(I) = 0D0
         Q (I) = 0D0
         QD(I) = 0D0
         end do
      if(k.eq.0) then
        CALL DUGOLM1(X,M,Nf,P,PD)
      IF(X.LT.1.D0) THEN
        A1=DSQRT(1.D0-X*X)
      ELSE
        A1=DSQRT(X*X-1.D0)
      ENDIF
      DO I = 1, Nf
        P(I)=P(I)*A1**M
        PD(I)=-PD(I)*A1**(M-2)
        IF(X.GT.1.D0) PD(I)=-PD(I)
      ENDDO
      end if
       ier1 = 0
      if(k.eq.1) CALL DLEGF1(P(1),D,X,D,d,M,M,K,0,IER1)
      CALL DLEGF2(Q(1),D,X,D,-M,M,K,0,IER2)

c        print *,m,q(1), x,d
c        write(7,*) m,x, q(1)
c        close(7)
c        pause

       IF((IER1+IER2).NE.0) WRITE(7,2) M,IER1,IER2
c       IF(IER2.NE.0) WRITE(7,2) M,IER2
    2 FORMAT(1X,'funlegnn  M,I,IER1,2,3,4=',6I10)
       ier3 = 0
       DO 3 I=2,NF
      if(k.eq.1) then
       if(i.gt.2) d=P(I-2)
      CALL DLEGF1(P(I),PD(I-1),X,P(I-1),d,M+I-1,M,K,1,IER3)
      end if
      CALL DLEGF2(Q(I),QD(I-1),X,Q(I-1),-M+I-1,M,K,1,IER4)

c        print *,m, i,q(i), q(i-1)
c        pause

c       IF((IER3+IER4).NE.0) WRITE(7,2) M,I,IER1,IER2,IER3,IER4
       IF(IER4.NE.0) WRITE(7,2) M,I,IER4
    3 CONTINUE
  900 RETURN
      END

c************************************************
C fact2
C
c      real FUNCTION fact2*8(n,k)
      real*8 FUNCTION fact2(n,k)
      REAL*8 f
      f = 1d0
       do i = k + 1, n
        f = f * i
       end do
      fact2 = f
      RETURN
      END

c************************************************
C fact3
C
      real*8 FUNCTION fact3 (n,k,f1)
      REAL*8 f, f1
      f = f1
       do i = k + 1, n
        f = f * i
       end do
      fact3 = f
      RETURN
      END

c************************************************
C ss1
C
      complex*16 FUNCTION ss1 (n)
      complex*16 s, s1, ss0
      REAL*8 aksi, ak
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      Ss0=(1D0,0D0)
      n1 = mod(iabs(n),4)
      go to (1, 2, 3, 4) n1+1
    1  ss1 = ss0
       RETURN
    2  ss1 = s1
       if(n.lt.0) ss1 = -ss1
       RETURN
    3  ss1 = -ss0
       RETURN
    4  ss1 = -s1
       if(n.lt.0) ss1 = -ss1
       RETURN
      END

c******************************************
      SUBROUTINE DUGOLM1(X,M,LMAX,PI,TAU)
      REAL*8 X,PI(LMAX),TAU(LMAX),A
      M2=2*M
      A=1.D0
      DO I=3,M2-1
       IF(MOD(I,2).EQ.1) A=A*I
      ENDDO
      PI(1)=A
      PI(2)=A*X*(M2+1.D0)
      TAU(1)=A*M*X
      TAU(2)=A*(M2+1.D0)*((M+1.D0)*X*X-1.D0)
      DO L1=3,LMAX
        PI(L1)=((2.D0*L1+M2-3.D0)*X*PI(L1-1)-(L1+M2-2.D0)*PI(L1-2))
     *          /(L1-1.D0)
      TAU(L1)=(L1+M-1.D0)*X*PI(L1)-(L1+M2-1.D0)*PI(L1-1)
      ENDDO
      RETURN
      END

c************************************************
C BESSJJ0
c
      SUBROUTINE BESSJJ0(A,NUM,BESJ)
c+&      parameter (nterm=200)
      IMPLICIT COMPLEX*16 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1)
   11 continue
c+&
c      BESJ(NUM+1)=0D0
          do i =1,num+1
            BESJ(i)=(0D0,0D0)
          end do
      BESJ(NUM)=(1D-300,1d-300)
      N=2*NUM+1
      NUM1=NUM-1
      DO I = 1, NUM1
       N=N-2
       I1=NUM-I
       BESj(I1) = n * a * BESj(I1+1) - BESj(I1+2)
         if(cdabs(besJ(i1)).gt.1d300) then
           num = i - 1
           go to 11
         end if
      end do
      N=2*(NUM/2)
      B=1.2533141373155002D0*CDSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
   13 BESJ(I)=C*BESJ(I)
      RETURN
      END

c************************************************
C CESSEL0
c
      SUBROUTINE CESSEL0(A,NUM,BESJ,BESY,iiy)
      IMPLICIT COMPLEX*16 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1),BESY(NUM+1),iiy(num+1)
      factor=1d250
   11 continue
          do i =1,num+1
            iiy(i) = 0
c+&
            BESJ(i)=(0D0,0D0)
            BESY(i)=(0D0,0D0)
          end do
c      BESJ(NUM+1)=(0D0,0D0)
      BESJ(NUM)=(1D-300,1D-300)
      N=2*NUM+1
      NUM1=NUM-1
      DO I=1,NUM1
       N=N-2
       I1=NUM-I
       BESj(I1) = n * a * BESj(I1+1) - BESj(I1+2)
        if(cdabs(besJ(i1)).gt.1d300) then
           num = i - 1
           go to 11
        end if
      end do
      B1=A*BESJ(1)-BESJ(2)
      B2=-A*B1-BESJ(1)
      N=2*(NUM/2)
      B=1.2533141373155002D0*CDSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
      BESJ(I)=C*BESJ(I)
   13 continue
c Y
      BESY(1)=-C*B1
      BESY(2)=C*B2
      DO 14 I=3,NUM
      I2=I-2
         if(iiy(i-1).eq.iiy(i2)) then
          BESY(I) = (2.0D0*I2+1.0D0) * a * BESY(I-1) - BESY(I2)
         else
          BESY(I) = (2.0D0*I2+1.0D0) * a * BESY(I-1) - BESY(I2)/factor
         end if
         iiy(i) = iiy(i-1)
         if(cdabs(besY(i)).gt.1d300) then
           besY(i) = besY(i) / factor
           iiy(i) = iiy(i-1) + 1
         end if
   14 continue
      RETURN
      END

c************************************************
C cdrkf45
c
      SUBROUTINE cDRKF45(F,NEQN,Y,T,TOUT,
     1   RELERR,ABSERR,IFLAG,WORK,IWORK)
      INTEGER IWORK(5),NEQN
c      REAL*8 Y(NEQN),T,TOUT,RELERR,ABSERR,
c     1   WORK(27)
      REAL*8 T,TOUT,RELERR,ABSERR
      complex*16 Y(NEQN), WORK(27)
      EXTERNAL F
      K1M = NEQN + 1
      K1 = K1M +1
      K2 = K1 + NEQN
      K3 = K2 + NEQN
      K4 = K3 + NEQN
      K5 = K4 + NEQN
      K6=K5+NEQN
      CALL cDRKFS (F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK(1),
     1   WORK(K1M),WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),
     2   WORK(K6),WORK(K6+1),IWORK(1),IWORK(2),IWORK(3),
     3   IWORK(4),IWORK(5))
      RETURN
      END

c************************************************
C cdrkfs
c
      SUBROUTINE cDRKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,
     1   IFLAG,YP,H,F1,F2,F3,F4,F5,SAVRE,SAVAE,
     2   NFE,KOP,INIT,JFLAG,KFLAG)
      LOGICAL  HFAILD, OUTPUT
      complex*16 Y(NEQN),YP(NEQN),F1(NEQN),
     *           F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN)
      REAL*8 T,TOUT,RELERR,ABSERR,H
c      REAL*8 Y(NEQN),T,TOUT,RELERR,ABSERR,H,YP(NEQN),F1(NEQN)
     1   ,SAVRE,SAVAE
      EXTERNAL F
      REAL*8  A,AE,DT,EE,EEOET,ESTTOL,ET,HMIN,REMIN,
     1   RER,S,SCALE,TOL,TOLN,U26,EPSP1,EPS,YPK
c   * added (Draine!)
      SAVE EPS,MAXNFE,REMIN,U26
      DATA EPS/1.D0/
c   *
      DATA REMIN /1.0D-16/
      DATA MAXNFE /500000/
c  *
      IF(EPS.GE.1.D0)THEN
    5    EPS=EPS/2.0D0
         EPSP1=EPS+1.0D0
         IF(EPSP1.GT.1.0D0)GOTO 5
         U26=26.0D0*EPS
      ENDIF
c  **************************************** (Draine!)
      IF(NEQN.LT.1) GO TO 10
      IF((RELERR.LT.0D0).OR.(ABSERR.LT.0D0)) GO TO 10
c        print *, iflag
c        pause
      MFLAG = IABS(IFLAG)
      IF((MFLAG.EQ.0).OR.(MFLAG.GT.8))  GO TO 10
      IF(MFLAG.NE.1) GO TO 20
c      EPS = 1.0D0
c    5 EPS = EPS/2.0D0
c      EPSP1 = EPS + 1.0D0
c      IF(EPSP1.GT.1.0D0) GO TO 5
c      U26 = 26.0D0 * EPS
      GO TO 50
   10 IFLAG =8
      RETURN
   20 IF((T.EQ.TOUT).AND.(KFLAG.NE.3)) GO TO 10
      IF(MFLAG.NE.2) GO TO 25
      IF((KFLAG.EQ.3).OR.(INIT.EQ.0)) GO TO 45
      IF(KFLAG.EQ.4) GO TO 40
      IF((KFLAG.EQ.5).AND.(ABSERR.EQ.0.0D0)) GO TO 30
      IF((KFLAG.EQ.6).AND.(RELERR.LE.SAVRE).AND.
     1   (ABSERR.LE.SAVAE)) GO TO 30
      GO TO 50
   25 IF(IFLAG.EQ.3) GO TO 45
      IF(IFLAG.EQ.4) GO TO 40
      IF((IFLAG.EQ.5).AND.(ABSERR.GT.0.0D0))  GO TO 45
   30 STOP
   40 NFE = 0
      IF(MFLAG.EQ.2) GO TO 50
   45 IFLAG = JFLAG
      IF(KFLAG.EQ.3)  MFLAG = IABS(IFLAG)
   50 JFLAG = IFLAG
      KFLAG = 0
      SAVRE = RELERR
      SAVAE = ABSERR
      RER = 2.0D0*EPS + REMIN
      IF(RELERR.GE.RER)  GO TO 55
      RELERR = RER
      IFLAG = 3
      KFLAG = 3
      RETURN
   55 DT = TOUT -T
      IF(MFLAG.EQ.1)  GO TO 60
      IF(INIT.EQ.0)  GO TO 65
      GO TO 80
   60 INIT = 0
      KOP = 0
      A = T
      CALL F(A,Y,YP)
      NFE = 1
      IF(T.NE.TOUT) GO TO 65
      IFLAG = 2
      RETURN
   65 INIT = 1
      H = DABS(DT)
      TOLN = 0.0D0
      DO 70 K = 1,NEQN
c         TOL=RELERR*DABS(Y(K)) + ABSERR
         TOL=RELERR*cDABS(Y(K)) + ABSERR
          IF(TOL.LE.0.0D0)  GO TO 70
         TOLN = TOL
c         YPK = DABS(YP(K))
         YPK = cDABS(YP(K))
         IF(YPK *H**5.GT.TOL)  H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF(TOLN.LE.0.0D0)  H = 0.0D0
      H=1D-4
      H= DMAX1(H, U26*DMAX1(DABS(T),DABS(DT)))
      JFLAG = ISIGN(2,IFLAG)
   80 H =DSIGN(H,DT)
      IF(DABS(H).GE.2.0D0*DABS(DT))  KOP = KOP+1
      IF(KOP.NE.100) GO TO 85
      KOP = 0
      IFLAG = 7
      RETURN
   85 IF(DABS(DT).GT.U26*DABS(T))  GO TO 95
      DO 90  K = 1,NEQN
   90 Y(K) = Y(K) + DT*YP(K)
      A = TOUT
      CALL F(A,Y,YP)
      NFE = NFE+ 1
      GO TO  300
   95 OUTPUT=.FALSE.
      SCALE= 2D0/RELERR
      AE = SCALE * ABSERR
  100 HFAILD = .FALSE.
      HMIN = U26 * DABS(T)
      DT = TOUT -T
      IF(DABS(DT).GE.2.0D0*DABS(H))  GO TO 200
      IF(DABS(DT).GT.DABS(H))  GO TO 150
      OUTPUT= .TRUE.
      H = DT
      GO TO  200
  150 H = 0.5D0 * DT
  200 IF(NFE.LE.MAXNFE)  GO TO 220
      IFLAG = 4
      KFLAG = 4
      RETURN
  220 CALL cFEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,F1)
      NFE = NFE + 5
      EEOET=0D0
      DO 250  K = 1,NEQN
c         ET = DABS(Y(K))+ DABS(F1(K)) + AE
         ET = cDABS(Y(K))+ cDABS(F1(K)) + AE
         IF(ET.GT.0.0D0)  GO TO 240
               IFLAG =5
               print *, iflag,neqn, k, ae, Y(K), F1(K), et
               pause
               RETURN
c  240    EE = DABS((-2.09D3 * YP(K) + (2.197D4 * F3(K) -
c     1       1.5048D4 * F4(K))) + (2.2528D4 *F2(K) -2.736D4 *F5(K)))
  240    EE = cDABS((-2.09D3 * YP(K) + (2.197D4 * F3(K) -
     1       1.5048D4 * F4(K))) + (2.2528D4 *F2(K) -2.736D4 *F5(K)))
  250    EEOET = DMAX1(EEOET,EE/ET)
      ESTTOL = DABS(H)*EEOET*SCALE/7.524D5
      IF (ESTTOL.LE.1.0D0)  GO TO 260
      HFAILD = .TRUE.
      OUTPUT = .FALSE.
      S = 0.1D0
      IF(ESTTOL.LT.5.9049D4)  S = 0.9D0/ESTTOL**0.2D0
      H = S * H
      IF(DABS(H).GT.HMIN)  GO TO  200
      IFLAG = 6
      KFLAG = 6
      RETURN
  260 T=T+H
      DO 270  K=1,NEQN
  270    Y(K) =F1(K)
      A = T
      CALL F(A,Y,YP)
      NFE = NFE + 1
      S = 5.0D0
      IF(ESTTOL.GT.1.889568D-4) S = 0.9D0/ESTTOL**0.2D0
      IF (HFAILD)  S = DMIN1(S,1.0D0)
      H= DSIGN( DMAX1(S*DABS(H),HMIN),H)
      IF(OUTPUT)  GO TO 300
      IF(IFLAG.GT.0) GO TO 100
      IFLAG = -2
      RETURN
  300 T = TOUT
      IFLAG = 2
      RETURN
      END

c************************************************
C cfehl
c
      SUBROUTINE cFEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
c      REAL*8  Y(NEQN),T,H,YP(NEQN),F1(NEQN),F2(NEQN),
c     1   F3(NEQN),F4(NEQN),F5(NEQN),S(NEQN),CH
      REAL*8  T,H,CH
      complex*16  Y(NEQN),YP(NEQN),F1(NEQN),S(NEQN),
     *            F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN)
      CH = H/4.0D0
      DO 221  K = 1,NEQN
  221    F5(K) = Y(K) + CH * YP(K)
      CALL F(T+CH,F5,F1)
      CH = 3.0D0 * H/3.2D1
      DO 222  K = 1,NEQN
  222    F5(K) = Y(K) + CH*(YP(K)+3.0D0*F1(K))
      CALL F(T+3.0D0*H/8.0D0,F5,F2)
      CH = H/2.197D3
      DO 223  K =1,NEQN
  223    F5(K) = Y(K)+CH*(1.932D3*YP(K)+(7.296D3*F2(K)-7.2D3*F1(K)))
      CALL F(T+1.2D1*H/1.3D1,F5,F3)
      CH = H/4.104D3
      DO 224  K= 1,NEQN
  224    F5(K)=Y(K)+CH*((8.341D3*YP(K)-8.45D2*F3(K))
     1         +(2.944D4*F2(K) - 3.2832D4*F1(K)))
      CALL F(T+H,F5,F4)
      CH = H/2.052D4
      DO 225 K= 1,NEQN
  225    F1(K)=Y(K)+CH*((-6.08D3*YP(K)+(9.295D3 *F3(K)-
     1         5.643D3*F4(K)))+(4.104D4*F1(K)-2.8352D4 *F2(K)))
      CALL F(T+H/2.0D0,F1,F5)
      CH = H/7.61805D6
      DO 230  K= 1,NEQN
  230    S(K)=Y(K)+CH*((9.0288D5*YP(K)+(3.855735D6 *F3(K) -
     1        1.371249D6*F4(K)))+(3.953664D6*F2(K)+2.7702D5*F5(K)))
      RETURN
      end

c************************************************
C cDRSF212
C
      SUBROUTINE cDRSF212(R,DR,ar,adr,i12,m,N,c,X,rl,AD,IA,BD,IB,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,S1,S,AD(IA),RL,cc,rll
      COMPLEX*16 WORK(27),Y(4),BD(IB), ar, adr, r, dr
      DIMENSION IWORK(5)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /F2/ RLl,Cc,Mm,i12a
      EXTERNAL  cDRF222
   11 FORMAT(1X,I10,2D10.3,4D20.8)
   31 FORMAT(17H TOLERANCES RESET,2D12.3)
   41 FORMAT(11H MANY STEPS)
   71 FORMAT(12H MUCH OUTPUT)
   81 FORMAT(14H IMPROPER CALL)

          mm = m
          cc = c
          rll = rl
          i12a = i12
c+&
      R=s
      Y(1) = s
      DR=s
      Y(2)=s
      aR=s
      Y(3)=s
      aDR=s
      Y(4)=s

      NEQN=2
      if(i12.eq.2) NEQN=4
       T=0D0
c       ELERR=1D-12
       ELERR=1D-15
      ABSERR=0D0
      TFINAL=X
      IFLAG=1
      TOUT=T
      CALL cDRF2K02(M,N,C,AD,IA,BD,IB,i12,Y)
c                             print *,c,neqn,i12
c        if(cdabs(c).lt.80d0)  print *,m,n,y(1),y(2), y(3), y(4)
c         if(n.ge.170)  print *,m,n,c,y(1),y(2), y(3), y(4)
c        if(n.ge.210)  print *,m,n,c,y(1),y(2), y(3), y(4)
c        if(m.eq.61.and.cdabs(c).lt.80d0)  print *,m,n,y(4)
   10 CALL cDRKF45(cDRF222,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG,WORK,IWORK)

c        if(m.eq.61.and.iflag.eq.5) then
c        if(m.eq.61.and.cdabs(c).lt.80d0) then
c         print *, m,n,iflag, ier
c         print *,  ELERR,ABSERR
c        print *,c
c         print *,y(1),y(2),y(3),y(4)
c         pause
c        end if

      if(i12.eq.0) then
      R=Y(1)
      DR=Y(2)
c                       print *,r,dr
      end if

      if(i12.eq.1) then
      aR=Y(1)
      aDR=Y(2)
      end if

      if(i12.eq.2) then
      R=Y(1)
      DR=Y(2)
      aR=Y(3)
      aDR=Y(4)
      end if

      GO TO (80,20,30,40,50,60,70,80),IFLAG
   20 TOUT=X
      IF(T.LT.TFINAL) GO TO 10
      IER=0
      GO TO  111
   30 WRITE(7,31)  ELERR,ABSERR
      IER=3
      GO TO 10
   40 WRITE(7,41)
      IER = 4
      GO TO 10
   50 ABSERR=1D-9
c       ELERR=10D0 * ELERR
      WRITE(7,31)  ELERR,ABSERR
      IER=5
c+&
c      iflag=-1
      GO TO 10
   60   ELERR=10D0 * ELERR
      WRITE(7,31)  ELERR,ABSERR
      IFLAG=2
      IER=6
c        print *, iflag, ier
c        pause
      GO TO 10
   70 WRITE(7,71)
      IFLAG=2
      IER=7
      GO TO 10
   80 WRITE(7,81)
      IER=8
      GO TO 1000
  111 CONTINUE
 1000 RETURN
c 1000 continue
c      print *, m,n,c,x
c      print *, iflag, ier
c      pause
c      RETURN
      END

c************************************************
C cDRF222
c
      SUBROUTINE cDRF222(X,Y,YP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C, RL, Y(4),YP(4), ar
      COMMON /F2/ RL,C,M,i12
       AKSI=X**2+1D0
       am = m * m
       ar = RL - (C*X)**2 - aM / AKSI

c        do k = 1,4 
c         YP(k)=(0d0,0d0)
c        end do 
         
      if(i12.eq.0.or.i12.eq.1) then
       YP(1)=Y(2)
       YP(2)=(-2D0*X*Y(2)+ aR * Y(1))/AKSI
c       YP(3)=(0d0,0d0)
c       YP(4)=(0d0,0d0)
      end if

      if(i12.eq.2) then
       YP(1)=Y(2)
       YP(2)=(-2D0*X*Y(2) + aR * Y(1))/AKSI
       YP(3)=Y(4)
       YP(4)=(-2D0*X*Y(4) + aR * Y(3))/AKSI
      end if

c       if(cdabs(yp(4)).gt.1d300.or.cdabs(y(4)).gt.1d300) then 
c          print *,m,y(4)
c          print *,m,yp(4)
c          pause
c        end if  

c       if(m.eq.61.and.x.gt.0.1005d0.and.cdabs(c).lt.80d0) then 
c       if(x.gt.0.1005d0.and.cdabs(c).lt.80d0) then 
c        print *,m,i12,c,x,rl
c        print *,y(1),y(2),y(3),y(4)
c        print *,yp(1),yp(2),yp(3),yp(4)
c        pause
c       end if
        
      RETURN
      END

c************************************************
C cDRF2K02
C
      SUBROUTINE cDRF2K02(M,N,C,AD,IA,BD,IB,i12,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AD,C,R1,S1,r,s,s2,t,r2,ss,ss2, bd, y
      DIMENSION AD(IA),BD(IB),FACT(300),Y(4)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /FACT/ FACT
      COMMON /eps3/ eps3
      COMMON /pi/ pi
       eps4 = 1d-40
c       eps4 = 1d-20
      M2=2*M
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 3
    2 JK=0
    3 CONTINUE
      S2 = s
      DO 4 I=1,IA
      J=2*(I-1)+JK
c      F=1D0
      F=ad(i)
      DO 5 JJ=1,M2
    5 F=F*(J+JJ)
c      T=AD(I)*F
      T=F
c         print *,i,s2,t
      S2=S2+T
      IF(cDABS(T/s2).LT.eps4) GO TO 6
    4 CONTINUE
    6 CONTINUE
      MF=N+M+1+JK
      F=FACT(MF)
                                
c+&
c      IF(MF.LE.169) F=FACT(MF)
c        IF(MF.ge.168) print *, m,n,mf,f,s2
            f=f/s2
c            IF(MF.ge.170) f=f*1d300
c      IF(MF.ge.170) stop 1706
c      IF(MF.ge.168) then
c        print *, m,n,mf,f,s2
c        pause
c      end if
      
c       print *, s2
c       print *, mf,f

c+&      if(i12.eq.1.or.i12.eq.2) then
c----
       jmn = mod(N-m-jk,4)
       if(jmn.eq.0) ss = 1d0
       if(jmn.eq.1) ss = s1
       if(jmn.eq.2) ss = -1d0
       if(jmn.eq.3) ss = -s1
c+&
       if(i12.eq.1.or.i12.eq.2) then
c       jmn2 = dabs(mod(m-n-jk,4))
       jmn2 = iabs(mod(m-n-jk,4))
       if(jmn2.eq.0) ss2 = 1d0
       if(jmn2.eq.1) ss2 = s1
       if(jmn2.eq.2) ss2 = -1d0
       if(jmn2.eq.3) ss2 = -s1
c-----
      if(mf.le.169)
     * R=Ss*(M2-1D0)*FACT(M+1)*PI/2D0**(2*N-M+1)*C**(M-1-JK)/
     *  FACT(M2+1)*
     *  (F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1))**2
     *  *s2/BD(M+1)
           if(mf.gt.169) then
       R=Ss*(M2-1D0)*FACT(M+1)*PI/2D0**(2*N-M+1)*C**(M-1-JK)/
     *  FACT(M2+1)*
     *  (F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1)*1d300)
     *  *s2/BD(M+1)
c          print *,r
        r=r*(F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1)*1d300)
       end if
c-----------------------------------------------------------------
      R1=Ss2*S2/2D0**M/C**(M+1+JK)/AD(1)/FACT(M+1)
c        if(m.eq.64.and.n.gt.207)  write(*,*) n,mf,
c     *  FACT(M2+1),BD(M+1),s2, r
      end if

      if(i12.eq.0.or.i12.eq.2) then
c          ss = 1d0
          R2=Ss*2D0**M*C**(M+JK)*AD(1)*FACT(M+1)/s2
c          write(7,*) c,n
c           if(n.eq.210)    write(*,*) s2,r2
c           if(m.eq.64)    write(*,*) c,n, s2,r2
c          pause
c          print *,AD(1),FACT(M+1),s2
      end if

      IF(JK.EQ.1) GO TO 10

      if(i12.eq.0) then
      Y(1)=R2/(m2+1d0)
      Y(2)=s
      end if

      if(i12.eq.1) then
      Y(1)=R
      Y(2)=R1*(M2+1D0)
      end if

      if(i12.eq.2) then
      Y(1)=R2/(m2+1d0)
      Y(2)=s
      Y(3)=R
      Y(4)=R1*(M2+1D0)
      end if

      GO TO 900

   10 CONTINUE

      if(i12.eq.0) then
      Y(1)=s
      Y(2)=R2/(M2+3D0)
      end if

      if(i12.eq.1) then
      Y(1)=R1*(M2+3D0)
      Y(2)=R*(M2-3D0)
      end if

      if(i12.eq.2) then
      Y(1)=s
      Y(2)=R2/(M2+3D0)
      Y(3)=R1*(M2+3D0)
      Y(4)=R*(M2-3D0)
c        if(cdabs(Y(4)).gt.1d260) then
c        if(n.ge.168) then
c           print *,m,n,jk,c,r,y(4)
c           print *,ss, fact(m+1), C**(M-1-JK)
c           print *,fact(m2+1), bd(m+1)
c           print *,f,s2
c           print *, mf, fact(mf)
c                 R=Ss*(M2-1D0)*FACT(M+1)*PI/2D0**(2*N-M+1)*C**(M-1-JK)/
c     *  FACT(M2+1)/BD(M+1)/S2*
c         print *, (F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1))
c           pause
c        end if
      end if

  900 CONTINUE

c      write(*,*) n,jk,c
c      write(*,22) y(1)
c      write(*,22) y(2)
c      write(7,22) y(1)
c      write(7,22) y(2)
c      write(*,22) y(3)
c      write(*,22) y(4)
c      write(7,22) y(3)
c      write(7,22) y(4)
c      pause
  22  format(1x,2d24.15)

      RETURN
      END

c************************************************
C cDRSF20
C
      SUBROUTINE cDRSF20(R,DR,X,R1,R2,m,c,rl,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,S1,S,RRA,RRB,R1,R2,R3,RL,cc,rll
      COMPLEX*16 r, dr, WORK(27),Y(4)
      DIMENSION IWORK(5)
      COMMON /K1/ S, S1, AKSI, AK, kkk, nk, nal
      COMMON /F2/ RLl, Cc, Mm, i12a
      EXTERNAL  cDRFP2
   11 FORMAT(1X,I10,2D10.3,4D20.8)
   31 FORMAT(17H TOLERANCES RESET,2D12.3)
   41 FORMAT(11H MANY STEPS)
   71 FORMAT(12H MUCH OUTPUT)
   81 FORMAT(14H IMPROPER CALL)
          mm = m
          cc = c
          rll = rl
      NEQN = 4
      PL = dreal(RL)
      AC = dreal(C)


c      print *, rL, C
c      print *, 1D0+PL/AC**2, 1D0-PL/AC**2
c      pause

      TKSI=DSQRT((1D0+PL/AC**2)/2D0+DSQRT((1D0-PL/AC**2)**2/4D0+
     *     (M**2-1D0)/AC**2))
      IF(M.EQ.1) TKSI=1D0
       T=TKSI+X/AC*10D0
*       ELERR=1D-5
      ELERR=1D-12
      ABSERR=0D0
      TFINAL=X
      IFLAG=1
      TOUT=T
      CALL cDRFP0(M,C,RL,T,Y)
   10 CALL cDRKF45(cDRFP2,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG,WORK,IWORK)
c      RRA=DCMPLX(Y(1),Y(2))
c      RRB=DCMPLX(Y(3),Y(4))
      RRA = Y(1) + Y(2)
      RRB = Y(3) + Y(4)
      R3=-S1/(C*AKSI*((RRB/RRA+X/AKSI)*R1+R2))
       R=-S1*R3
       DR=-(X/AKSI+RRB/RRA)*R+R1*S1*RRB/RRA
      GO TO (80,20,30,40,50,60,70,80),IFLAG
   20 TOUT=X
      IF(T .GT.TFINAL) GO TO 10
      IER=0
      GO TO  111
   30 WRITE(7,31)  ELERR,ABSERR
      IER=3
      GO TO 10
   40 WRITE(7,41)
      IER = 4
      GO TO 10
   50 ABSERR=1D-9
      WRITE(7,31)  ELERR,ABSERR
      IER=5
      GO TO 10
   60  ELERR=10D0* ELERR
      WRITE(7,31)  ELERR,ABSERR
      IFLAG=2
      IER=6
      GO TO 10
   70 WRITE(7,71)
      IFLAG=2
      IER=7
      GO TO 10
   80 WRITE(7,81)
      IER=8
      GO TO 1000
  111 CONTINUE
 1000 RETURN
      END

c************************************************
C cDRFP0
C
      SUBROUTINE cDRFP0(M,C,RL,KSI,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 RL,C, y, c2, gksi, al, an
      REAL*8 KSI
      DIMENSION Y(4)
      aM2 = M * M
      C2=C**2
      AKSI=KSI**2-1D0
      GKSI=C2-(RL-C2)/AKSI+(1D0-aM2)/AKSI**2
      AL=KSI*(2D0*(1D0-aM2)-(RL-C2)*AKSI)/(2D0*GKSI*AKSI**3)
      AN=cDSQRT(1D0+cdabs(GKSI)+AL**2)
      Y(1)=1D0/AN
      Y(2)=0D0
      Y(3)=-AL/AN
      Y(4)=-DSQRT(cdabs(GKSI))/AN
      RETURN
      END

c************************************************
C cDRFP2
c
      SUBROUTINE cDRFP2(X,Y,YP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,RL, Y(4),YP(4), c2, gksi, akk
      COMMON /F2/ RL,C,M, i1
      AKSI=X**2-1D0
      C2=C**2
      GKSI=C2-(RL-C2)/AKSI+(1D0-M**2)/AKSI**2
      AKK=(Y(1)*Y(3)+Y(2)*Y(4)*(1D0-GKSI)/(Y(1)**2+Y(2)**2+Y(3)**2+
     *     Y(4)**2))
       YP(1)=AKK*Y(1)-Y(3)
       YP(2)=AKK*Y(2)-Y(4)
       YP(3)=AKK*Y(3)+GKSI*Y(1)
       YP(4)=AKK*Y(4)+GKSI*Y(2)
      RETURN
      END

c************************************************
c CDRG1bn
c
      SUBROUTINE CDRG1bn(R,DR,AR,ADR,I12,M,n,C,x,AD,IA,BD,IB,rrl,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 C, dr, drr, drr2, dr2, AD(IA), drr3
c      COMPLEX*16 ar, adr
      DIMENSION SG(1000), sa(1000), sb(1000), BD(IB+1),
     *          sam(100)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /FACT/ FACT(300)
      COMMON /eps3/ eps3
      ier=0
c      if(i12.eq.0) go to 900
c! lambda according Slavyanov
      rl = rrl - c**2
c+&      eps4 = 1d-20
      eps4 = 1d-100
      aM = 0.5d0 * M
      if(k.eq.0) xx0 = ((x - 1d0) / (x + 1d0)) ** am
      if(k.eq.0) xx =  (x - 1d0) / (x + 1d0)
      rp = - S1 * c
      re = cdexp(-rp * (x - 1d0))
      rp2 = 2d0 * rp
      rrp2l = (rp2 + 1d0) * (m + 1d0) + rl

c coefficients
      M2=2*M
      nm1=n-m
      nm2=n+m
      IF(MOD((N-M),2)) 1,2,1
    1  CONTINUE
      JK=1
      GO TO 36
    2  CONTINUE
      JK=0
   36 CONTINUE

c       if((nm1-jk)/2+1.gt.170) stop 1707
      ff3 = fact3(nm2+jk,(nm2+jk)/2,1d0/fact((nm1-jk)/2+1))
      rg0=(2d0*m+1d0+jk*2d0)/2d0**(nm2+1d0) * ff3 / AD(1)/C**(M+JK)
      rg0 = 1d0 / rg0

c 2===
      IF(I12.ne.0) then
      IF(MOD(m,2)) 33,34,33
   33  CONTINUE
      Jm=1
      GO TO 35
   34  CONTINUE
      Jm=0
   35 CONTINUE
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
      rf2 = AD(I) * fact2(m2+j,j)
      SIG=SIG+rf2
      if(cdabs(rf2/sig).lt.eps4) go to 8
    7  CONTINUE
    8  CONTINUE
      RK2=2d0**nm1/(m2-1d0)*fact(m2+1)/fact(m+1) / ff3
     *    * bd(m+1)*sig/C**(M-1-JK)
      if(jk.eq.1.and.k.eq.0) rk2=-rk2/(m2-3d0)
      if(k.eq.0) go to 339
      if(jk.eq.1.and.jm.eq.0) rk2=rk2/(m2-3d0)
      if(jk.eq.1.and.jm.eq.1) rk2=-rk2/(m2-3d0)*(-s1)
      if(jk.eq.0.and.jm.eq.0) rk2=rk2*s1
      if(k.eq.1) rk2=rk2/(-s1)**(M-1-JK)
  339 continue
      end if
c 2===
c------------
      Sssa0 = S
c. 1s
      DO I = 1, M
        j = - 2 * (M - I + 1) + JK
        if(j.ge.-m) go to 351
      sssa0 = sssa0 + BD(m-I+2) * (-1)**(-j-1) * fact(-j)
     *        * fact(m2+j+1)
      end do

351   continue

c. 2s
      rf1 = s
      DO I = 1, IA
      ir = 2* (i - 1) + jk
      rf3 = s
c=
         do j = 1, m
           ff2 = fact(m+1) / fact(m-j+1)
           rf2 = s
           do jj = 1, j+1
            jjj = jj - 1
            rf2 = rf2 + fact2(m,m-jjj)**2 / fact(jjj+1) /
     *            fact(j-jjj+1) * fact2(ir+m2-jjj,ir+jjj)
           end do
           rf2 = rf2 * (-1)**j / (2d0 * j)  * ff2
           rf3 = rf3 + (-1)**(j-1) / (2d0 * j)  * ff2**2 *
     *           fact2(ir+m2-j,ir+j) + rf2
         end do
c=
         irm = (ir + m - 1 + jk) / 2
         rf4 = s
         if(ir.eq.0) go to 295
         do j = 1, irm
           jj = j - 1
           if(2*jj.gt.ir-1) go to 295
           rf4 = rf4 + (2d0 * ir + m2 - 4d0 * jj - 1d0) / (ir + m - jj)
     *           / (2d0 * jj + 1d0) * fact2(ir+m2-2*jj-1,ir-2*jj-1)
         end do
 295   CONTINUE
      rf11 =  AD(I) * (rf3 - rf4)
      rf1 = rf1 + rf11
       if(cdabs(rf11/rf1).lt.eps4) go to 297
      end do
 297   CONTINUE
      Sssa0 = Sssa0 + rf1
      rs2 = rf1

c. 3s
      rf1 = S
      MR = M + 1
      DO I = MR, IB
        J = 2 * I  - JK
        rf11 = bD(I+1) * fact2(j-1,j-m2-1)
        rf1 = rf1 + rf11
c      if(cdabs(rf11).lt.1d-100) go to 298
      if(cdabs(rf11/rf1).lt.eps4) go to 298
      end do

  298 continue

      Sssa0 = Sssa0 + rf1

c 1--
      SG0 = (1d0,0d0)
      SG(1) = rrp2l / (m + 1d0) * sg0
      SG(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sg(1) -
     *        (m + 1d0) * sg0) / (2d0 * (m + 2d0))
      r = sg0 + sg(1) * xx + sg(2) * xx**2

c       print *,r
c       write(7,*)  r

c      r = (sg0 + sg(1) * xx + sg(2) * xx**2) * re
      drr = sg(1) * 1d0 +  sg(2) * 2d0 * xx
c      drr = (sg(1) * 1d0 +  sg(2) * 2d0 * xx) * re

c 2--
      if(i12.ne.0) then

      do i = 1, 100
         Sam(i) = s
      end do

      Sam(1) = - 2d0 / m / c / rg0
      if(-m+1.eq.0) go to 352
      Sam(2) = (2d0 * m * (rp2 + 1d0) - rrp2l) / (m - 1d0) * sam(1)
      if(-m+2.eq.0) go to 352
       do i = 1, m-2
         j = - m + i
         Sam(i+2) = ((2d0 * j * (j + rp2 + m + 1d0) + rrp2l)
     *              * Sam(i+1) -
     *           (j * (j + m) * sam(i))) / ((j + 1d0) * (j + m + 1d0))
       end do

  352 continue

c a0--
      s2ic = 2d0 * s1 * c
      sssa1 = s
       do i = 1, m
       j = i - 1
       rf = s
         do ii = 1, m - j
           rf = rf + s2ic**ii / fact(ii+1) * fact(m-j) / fact (ii)
     *          / fact(m-ii-j+1)
         end do
       sssa1 = sssa1 + sam(i) * rf
       end do

      Sa00 = 2d0 / fact(m+1) * sssa0 / rk2 - sssa1

c a0--
c!!!!!!!!          Sa0 = sa00

      Sb0 = (-2d0 * (rp2 + m) + rrp2l) * Sam(m) / m
      if(m.gt.1) sb0 = sb0 + (m - 1d0) * sam(m-1) / m

      sb(1) = rrp2l / (m + 1d0) * sb0
      sb(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sb(1) -
     *        (m + 1d0) * sb0) / (2d0 * (m + 2d0))
      Sa(1) = (rrp2l * sa00 - (m + 2d0) * sb(1)
     *        + 2d0 * (rp2 + m + 1d0) * sb0) / (m + 1d0)
      Sa(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sa(1)
     *        - (m + 1d0) * sa00 - ((m + 4d0) * sb(2)
     *        - 2d0 * (rp2 + m + 3d0) * sb(1)
     *        + (m + 2d0) * sb0)) / (2d0 * (m + 2d0))

      r2 = s
      drr2 = s
       do i = 1, m
         r2 = r2 + sam(m-i+1) / xx**i
         drr2 = drr2 - i * sam(m-i+1) / xx**(i+1)
       end do

      r2 = r2 + sa00 + sa(1) * xx + sa(2) * xx**2
      rr2 = sb0 + sb(1) * xx + sb(2) * xx**2
      drr2 = drr2 + sa(1) * 1d0 +
     *       2d0 * sa(2) * xx
      drr3 = sb(1) * 1d0 + 2d0 * sb(2) * xx
c      r2 = r2 + sa00 + sa(1) * xx + sa(2) * xx**2
c      r2 = r2 * re
c      rr2 = (sb0 + sb(1) * xx + sb(2) * xx**2) * re
c      drr2 = drr2 + sa(1) * 1d0 +
c     *       2d0 * sa(2) * xx
c      drr2 = drr2 * re
c      drr3 = (sb(1) * 1d0 + 2d0 * sb(2) * xx) * re
      end if

c      cosn = cdcos(c*(x-1d0))
c      csin = cdsin(c*(x-1d0))

c          write(7,*)  r

      DO  I = 2, 998
c 1--
      a = (i + 1d0) * (i + m + 1d0)
      rb = 2d0 * i * (i + m + 1d0 + rp2) + rrp2l
      g = i * (i + m)

      sg(i+1) = (rb * sg(i) - g * sg(i-1)) / a
      rs = sg(i+1) * xx**(i+1d0)
c      drs = dreal(rs)
c      airs = dimag(rs)
c      rsn = drs * cosn - airs * csin + s1 * (drs * csin + airs * cosn)

c      if(dreal(rs).ge.0d0) Rr01 = Rr01 + RS
c      if(dreal(rs).lt.0d0) Rr02 = Rr02 + RS
      R = R + RS

c          write(7,*) i, r, rsn
c          write(7,*) rs

c       R = R + RSn

      rd = sg(i+1) * (i + 1d0) * xx**i
c      drn = dreal(rd)
c      airn = dimag(rd)
c      rdn = drn * cosn - airn * csin + s1 * (drn * csin + airn * cosn)

      dRr = dRr + Rd
c      dRr = dRr + Rdn
c 2--
      rs2 = s
      rrs2 = s
      rd21 = s
      rd23 = s

      if(i12.ne.0) then

        if(cdabs(sb(i)).gt.1d300) return

      sb(i+1) = (rb * sb(i) - g * sb(i-1)) / a
      sa(i+1) = (rb * sa(i) - g * sa(i-1) - ((2d0 * i + m + 2d0) *
     *          sb(i+1) - 2d0 * (2d0 * i + rp2 + m + 1d0) * sb(i) +
     *          (2d0 * i + m) * sb(i-1))) / a
      rs2 = sa(i+1) * xx**(i+1d0)
      rrs2 = sb(i+1) * xx**(i+1d0)

c      drs2 = dreal(rs2)
c      airs2 = dimag(rs2)
c      rs2n = drs2 * cosn - airs2 * csin + s1 * (drs2 * csin
c     *       + airs2 * cosn)
c
c      drs22 = dreal(rrs2)
c      airs22 = dimag(rrs2)
c      rrs2n = drs22 * cosn - airs22 * csin + s1 * (drs22 * csin
c     *       + airs22 * cosn)

      R2 = R2 + RS2
c      R2 = R2 + RS2n
      Rr2 = Rr2 + RrS2
c      Rr2 = Rr2 + RrS2n

      rd21 = sa(i+1) * (i + 1d0) * xx**i
      rd23 = sb(i+1) * (i + 1d0) * xx**i
c      drd21 = dreal(rd21)
c      aird21 = dimag(rd21)
c      rd21n = drd21 * cosn - aird21 * csin + s1 * (drd21 * csin
c     *       + aird21 * cosn)
c      drd23 = dreal(rd23)
c      aird23 = dimag(rd23)
c      rd23n = drd23 * cosn - aird23 * csin + s1 * (drd23 * csin
c     *       + aird23 * cosn)

      dRr2 = dRr2 + rd21
      dRr3 = dRr3 + rd23
c      dRr2 = dRr2 + rd21n
c      dRr3 = dRr3 + rd23n

c      if(CDABS(RS2n/r2).LT.eps4.AND.CDABS(Rrs2n/rr2).LT.eps4) GO TO 15
      if(CDABS(RS2/r2).LT.eps4.AND.CDABS(Rrs2/rr2).LT.eps4) GO TO 15
      end if
   15 CONTINUE
c      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RDn/drr).LT.eps4) GO TO 14
      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RD/drr).LT.eps4) GO TO 14
      end do
   14 CONTINUE
c 1--

c       R = R + rr01 + rr02

c      print *,r

      R = R * xx0 / (x + 1d0) * re * rg0


c      R = R * xx0 / (x + 1d0) * rg0

c      print *,r
c      print *,xx0, x, re, rg0

c      aaR = R
c      R = aaR
      dr = r * (-rp  - 1d0/(x + 1d0) + m / aksi) +
     *     DRr * (2d0 / (x + 1d0)**3) * xx0 * re * rg0
c      dr = r * (-rp  - 1d0/(x + 1d0) + m / aksi) +
c     *     DRr * (2d0 / (x + 1d0)**3) * xx0 * rg0
c      aaR = dR
c      dR = aaR
c 2--
      if(i12.ne.0) then
      R2 = (rr2 * dlog(xx) + R2) * xx0 / (x + 1d0) * re
c      R2 = (rr2 * dlog(xx) + R2) * xx0 / (x + 1d0)
      dr2 = r2 * (-rp -1d0/(x + 1d0) + m / aksi) +
     *     (DRr2 + rr2 / xx + drr3 * dlog(xx)) *
     *     (2d0 / (x + 1d0)**3) * xx0 * re
c      dr2 = r2 * (-rp -1d0/(x + 1d0) + m / aksi) +
c     *     (DRr2 + rr2 / xx + drr3 * dlog(xx)) *
c     *     (2d0 / (x + 1d0)**3) * xx0
      AR = r2
      ADR = dr2
      end if

      if(dimag(c).eq.0d0) then
      r = dreal(r)
      dr = dreal(dr)
      end if

  900 RETURN
      END

c
c        SUBROUTINE SEGV(M,N,C,KD,CV,EG)
        SUBROUTINE SEGV1(M,jN,rC1,KD,CV,rlc1)
C
C       =========================================================
C       Purpose: Compute the characteristic values of spheroidal
C                wave functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for Prolate; KD=-1 for Oblate
C       Output:  CV --- Characteristic value for given m, n and c
C                EG(L) --- Characteristic value for mode m and n'
C                          ( L = n' - m + 1 )
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION B(1000),H(1000),D(3000),E(3000),F(3000),CV0(1000),
     &            A(3000),G(3000),EG(2000)
        complex*16 rc1, rlc1(2000)
c        DIMENSION B(100),H(100),D(300),E(300),F(300),CV0(100),
c     &            A(300),G(300),EG(200)
        c = rc1
        n = jn + m - 1

c            write(*,*) m,jn, n,kd, c
c            pause

        IF (C.LT.1.0D-10) THEN
           DO 5 I=1,N
5             EG(I)=(I+M)*(I+M-1.0D0)
           GO TO 70
        ENDIF                                           
        ICM=(N-M+2)/2
        NM=10+INT(0.5*(N-M)+C)
        CS=C*C*KD
        DO 60 L=0,1
           DO 10 I=1,NM
              IF (L.EQ.0) K=2*(I-1)
              IF (L.EQ.1) K=2*I-1
              DK0=M+K
              DK1=M+K+1
              DK2=2*(M+K)
              D2K=2*M+K
              A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
              D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &             *(DK2+3.0))*CS
10            G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
           DO 15 K=2,NM
              E(K)=DSQRT(A(K-1)*G(K))
15            F(K)=E(K)*E(K)
           F(1)=0.0D0
           E(1)=0.0D0
           XA=D(NM)+DABS(E(NM))
           XB=D(NM)-DABS(E(NM))
           NM1=NM-1
           DO 20 I=1,NM1
              T=DABS(E(I))+DABS(E(I+1))
              T1=D(I)+T
              IF (XA.LT.T1) XA=T1
              T1=D(I)-T
              IF (T1.LT.XB) XB=T1
20         CONTINUE
           DO 25 I=1,ICM
              B(I)=XA
25            H(I)=XB
           DO 55 K=1,ICM
              DO 30 K1=K,ICM
                 IF (B(K1).LT.B(K)) THEN
                    B(K)=B(K1)
                    GO TO 35
                 ENDIF
30            CONTINUE
c35            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
35            IF (K.NE.1) then
                 IF(H(K).LT.H(K-1)) H(K)=H(K-1)
              ENDIF
40            X1=(B(K)+H(K))/2.0D0
              CV0(K)=X1
c+&              IF (DABS((B(K)-H(K))/X1).LT.1.0D-14) GO TO 50
              IF (DABS((B(K)-H(K))/X1).LT.1.0D-15) GO TO 50
              J=0
              S=1.0D0
              DO 45 I=1,NM
                 IF (S.EQ.0.0D0) S=S+1.0D-30
                 T=F(I)/S
                 S=D(I)-T-X1
                 IF (S.LT.0.0D0) J=J+1
45            CONTINUE
              IF (J.LT.K) THEN
                 H(K)=X1
              ELSE
                 B(K)=X1
                 IF (J.GE.ICM) THEN
                    B(ICM)=X1
                 ELSE
                    IF (H(J+1).LT.X1) H(J+1)=X1
                    IF (X1.LT.B(J)) B(J)=X1
                 ENDIF
              ENDIF
              GO TO 40
50            CV0(K)=X1
              IF (L.EQ.0) EG(2*K-1)=CV0(K)
              IF (L.EQ.1) EG(2*K)=CV0(K)
55         CONTINUE
60      CONTINUE
70      CONTINUE
c70      CV=EG(N-M+1)
        CV=EG(N-M+1)
        rlc1(n-m+1) = CV

c            write(*,*) m,n, cv
c            write(*,*) n-m+1
c            write(*,*)  rlc1(n-m+1)
c            pause

        RETURN
        END

c
c      subroutine CDLAM_m(RLC2,M,NN,Cin,EPS,IEr)
      subroutine CDLAM_m(RLC2,M,NN,Cin,IEr)
c complex 'lambda' according Martin....
c vector of eigenvalues found from diagonalizing real symmetric tri-diag
c matrix(tql1).finds egv(m,m),egv(m+1,m)... as xlam(1),xlam(2)... out
c and in -side spheroid.accuracy is improved and d(r,mn) made by bowkmp
      parameter (nterms=330)
      implicit real*8(a-h,o-z)
c ** s1
      complex*16 cin, cdiag,csubd,cc,ss,rs1,rlc2(nterms)
      COMMON /K1/ sS, rS1, AKSI, AK, K, NK, nal
      dimension cdiag(nterms),csubd(nterms)
      integer r,s
      maxdim=nterms
c      maxdim=100
      msq=m**2
      m2=2*m
      sgn=ak
c      nl=kmax1
      nl=nn
      numax=nl-1
cm* now do inside
      cc=cin
c      print *, M,NN, cin
c      pause
c find eigenvalues of even parity first
      s=0d0
c determine matrix size and check for overflow and other errors
      n=numax+3*idint(dsqrt(cdabs(cc)))+2
      if((nl.le.0).or.(n.gt.maxdim)) goto 100
c set the diagonal elements.note that we write r=2*q+s
 1005 r=s
      do 1001 i=1,n
         xa=dfloat((m+r)*(m+r+1))
         cdiag(i)=(2.0d0*xa-dfloat(2*msq+1))*sgn*cc*cc
     $        /dfloat((m2+2*r-1)*(m2+2*r+3))+xa
 1001 r=r+2
c now set the subdiagonal elements
      r=s+2d0
      do 1002 i=2,n
         csubd(i)=dsqrt(dfloat((m2+r)*(m2+r-1))*dfloat(r*(r-1))
     $        /dfloat((m2+2*r-3)*(m2+2*r+1)))*cc*cc/dfloat(m2+2*r-1)
 1002 r=r+2
      ib=0
      call sf02(n,cdiag,csubd,ib)
cm* this goes to setup for and execution of nag routine f02ajf and sorts
      if(ib.ne.0) goto 102
c sf02 returns eigenvalues of same parity in ascending order,
cm* with respect to the size of the real part (see onaka, asano's letter)
cm? is this completely general?
      j=0
      ls=s+1
      do 1030 i=ls,nl,2
         rlc2(i)=cdiag(j+1)
c      write(6,*) i,rlc2(i)
c      print *, i,rlc2(i)
c      pause
 1030 j=j+1
      if(s.eq.1) return
c find eigenvalues of odd parity
      s=1
      goto 1005
cm*
cm* error returns
 100  if(nl.le.0) write(6,600)
      if (n.gt.maxdim) write(6,602) n
c+&
      if (n.gt.maxdim) n=maxdim
c      if (n.gt.maxdim) stop 100
 600  format(20h-nl negative or zero)
 602  format(' -allowable array size exceeded. increase dimension 
     $   statements in sphwf to at least ',i4)
      return
 102  write(6,604) ib
 604  format(20h-f02ajf error number,i4)
      ier=ib
      return
      end

c
      subroutine sf02(n,cdiag,csubd,ifail)
      parameter (nterms=330)
      implicit real*8(a-h,o-z)
      complex*16 cdiag(n),csubd(n),dcmplx
c      dimension ar(100,100),ai(100,100),rr(100),ri(100)
c      dimension intger(100),irank(100)
      dimension ar(nterms,nterms),ai(nterms,nterms),rr(nterms),
     *          ri(nterms),intger(nterms),irank(nterms)
      do 1 i=1,n
         do 1 j=1,n
            ar(i,j)=0.
            ai(i,j)=0.
 1    continue
cm* diagonal
      do 2 i=1,n
         ar(i,i)=dreal(cdiag(i))
         ai(i,i)=dimag(cdiag(i))
 2    continue
cm* subdiagonal, and superdiagonal since symmetric
      do 3 i=2,n
         arr=dreal(csubd(i))
         brr=dimag(csubd(i))
         j=i-1
         ar(i,j) = arr
         ar(j,i) = arr
         ai(i,j) = brr
         ai(j,i) = brr
 3    continue
      ifail = 1
      call f02ajf(ar,nterms,ai,nterms,n,rr,ri,intger,ifail)
      if(ifail.eq.0) go to 20
      write(6,*) ' f02ajf',ifail
      stop
cm* rank the real parts, then rearrange, using mark 12 nag routines
c
c 20   continue
c      ifail=0
c      call m01daf(rr,1,n,'a',irank,ifail)
ccm* rank the real parts rr using a mark 1 nag routine
  20   ifail=1
       call m01aaf(rr,1,n,irank,intger,ifail)
      if(ifail.eq.0) go to 21
c      write(6,*) ' m01daf ',ifail
      write(6,*) ' m01aaf ',ifail
      stop
cm* sort into complex number
c in sub. sf02
c!!! ordering - here! cdiag(irank(i))-YES,cdiag(i)-NO
 21   do 9 i=1,n
c         cdiag(irank(i)) =dcmplx(rr(i),ri(i))
c          print *, i,irank(i)
c          pause
         cdiag(i) =dcmplx(rr(i),ri(i))
 9    continue
c
c      if(.not.wanted) return
c      write(6,99994)
c99994 format(1x,42hcomplex eigenvalues - even then odd parity)
c      write(6,99995) (cdiag(i),i=1,n)
c      do i=1,n
c        write(6,*) i, cdiag(i)
c      end do
c99995 format(1x,3(2e12.4,2x))
      return
      end

c
      SUBROUTINE F02AJF(AR, IAR, AI, IAI, N, RR, RI, INTGER, LFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C
C     EIGENVALUES OF A COMPLEX*16 UNSYMMETRIC MATRIX
C     1ST AUGUST 1971
C
      INTEGER P01AAF, ISAVE, LFAIL, N, IAR, IAI, INTGER(N)
C$P 1
c * Draine!
      CHARACTER*8 SRNAME
      SAVE SRNAME
c      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION ACC, AR(IAR,N), AI(IAI,N),	RR(N), RI(N),
     * X02AAF
      DATA SRNAME /8H F02AJF /
      ISAVE = LFAIL
      LFAIL = 1
      ACC = X02AAF(1.0D0)
      CALL F01AMF(N, 1,	N, AR, IAR, AI,	IAI, INTGER)
      CALL F02ANF(N, ACC, AR, IAR, AI, IAI, RR,	RI, LFAIL)
      IF (LFAIL.NE.0) LFAIL = P01AAF(ISAVE,LFAIL,SRNAME)
      RETURN
      END

c
      SUBROUTINE M01AAF(A, MM, NN, IP, IST, IFAIL)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 10B REVISED.	IER-409	(JAN 1983).
C     DETACHED KEY SORT, TO SORT AN INDEX ARRAY	TO GIVE	THE
C     ASCENDING	ORDER OF A REAL	ARRAY.
      INTEGER P01AAF, M, MM, N,	NN, I, IFAIL, LIMIT, JS, L, J, M2,K,
     * IP(NN), IST(NN)
C$P 1
c * Draine!
      CHARACTER*8 SRNAME
      SAVE SRNAME
c      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION A(NN)
      DATA SRNAME /8H M01AAF /
      M	= MM
      N	= NN
      IF (MM) 440, 440,	20
   20 IF (NN) 420, 420,	40
   40 IF (NN-MM) 460, 400, 60
   60 DO 80 I=M,N
	 IP(I) = I
	 IST(I)	= I
   80 CONTINUE
      IFAIL = 0
  100 LIMIT = (N+M-1)/2
      JS = M - 2
      L	= M - 1
  120 L	= L + 1
      JS = JS +	2
      I	= IST(JS)
      J	= IST(JS+1)
  140 IF (A(I).GT.A(J))	GO TO 160
      M2 = I
      I	= IP(I)
      GO TO 180
  160 M2 = J
      J	= IP(J)
  180 K	= M2
  200 IF (A(I).GT.A(J))	GO TO 220
      IP(K) = I
      K	= I
      I	= IP(I)
      IF (IP(K)-K) 200,	260, 200
  220 IP(K) = J
      K	= J
      J	= IP(J)
      IF (IP(K)-K) 200,	240, 200
  240 J	= I
  260 IP(K) = J
      K	= J
      IF (IP(J)-J) 280,	300, 280
  280 J	= IP(J)
      GO TO 260
  300 IST(L) = M2
      IF (L.LT.LIMIT) GO TO 120
      IF (JS+1-N) 320, 340, 340
  320 JS = JS +	2
      I	= IST(L)
      J	= IST(JS)
      GO TO 140
  340 N	= LIMIT
      IF (LIMIT.GT.M) GO TO 100
      I	= IST(M)
      J	= 0
  360 IF (IP(I).EQ.I) GO TO 380
      J	= J + 1
      K	= IP(I)
      IP(I) = J
      I	= K
      GO TO 360
  380 IP(I) = J	+ 1
      GO TO 480
  400 IP(MM) = 1
      IFAIL = 0
      GO TO 480
  420 IFAIL = P01AAF(IFAIL,3,SRNAME)
      GO TO 480
  440 IFAIL = P01AAF(IFAIL,2,SRNAME)
      GO TO 480
  460 IFAIL = P01AAF(IFAIL,1,SRNAME)
  480 RETURN
      END

c
      SUBROUTINE F01AMZ(AR, IAR, AI, IAI, M, N,	BR, IBR, BI, IBI,CR, CI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C
C     COMPUTES	C = C +	 A*B  (COMPLEX)	WHERE
C     A	IS RECTANGULAR M BY N.
C     C	MUST BE	DISTINCT FROM B.
C     THE ELEMENTS OF B	MAY BE NON-CONSECUTIVE,	WITH OFFSETS IBR
C     AND IBI.
C
C     .. SCALAR	ARGUMENTS ..
      INTEGER IAR, IAI,	IBR, IBI, M, N
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION AR(IAR,N), AI(IAI,N), BR(IBR,N),	BI(IBI,N), CR(M)
     *,	CI(M)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I, J
C     ..
C
      DO 60 J=1,N
	 DO 40 I=1,M
	    CR(I) = (CR(I)+AR(I,J)*BR(1,J)) - AI(I,J)*BI(1,J)
	    CI(I) = (CI(I)+AR(I,J)*BI(1,J)) + AI(I,J)*BR(1,J)
   40	 CONTINUE
   60 CONTINUE
      RETURN
      END

c
      SUBROUTINE F02ANF(NN, ACC, HR, IHR, HI, IHI, WR, WI, LFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 9 REVISED. IER-326 (SEP 1981).
C
C     COMLR
C     FINDS THE	EIGENVALUES OF A COMPLEX UPPER-HESSENBERG MATRIX,
C     H, STORED	IN THE ARRAYS HR(N,N) AND HI(N,N), AND
C     STORES THE REAL PARTS IN THE ARRAY WR(N) AND THE COMPLEX
C     PARTS IN THE ARRAY WI(N).	ACC IS THE RELATIVE MACHINE
C     PRECISION. THE SUBROUTINE	FAILS IF ALL EIGENVALUES TAKE MORE
C     THAN 30*N	ITERATIONS.
C     1ST DECEMBER 1971
C
      INTEGER P01AAF, ISAVE, LFAIL, N, NN, ITS,	L, LL, N1, I, J, M,MM,
     * M1, I1, IHR, IHI, ITN
C$P 1
c * Draine!
      CHARACTER*8 SRNAME
      SAVE SRNAME
c      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION TR, TI, ACC, SR,	SI, XR,	XI, YR,	YI, ZR,	ZI, ZXR,
     * ZXI,HR(IHR,NN), HI(IHI,NN), WR(NN), WI(NN)
      DATA SRNAME /8H F02ANF /
      ISAVE = LFAIL
      LFAIL = 0
      N	= NN
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
   20 IF (N.EQ.0) GO TO	480
      ITS = 0
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
   40 L	= N + 1
      IF (N.EQ.1) GO TO	80
      DO 60 LL=2,N
	 L = L - 1
	 IF ((DABS(HR(L,L-1))+DABS(HI(L,L-1))).LE.(ACC*(DABS(HR(L-1,L-1)
     *	 )+DABS(HI(L-1,L-1))+DABS(HR(L,L))+DABS(HI(L,L))))) GO TO100
   60 CONTINUE
   80 L	= 1
  100 IF (L.EQ.N) GO TO	460
      IF (ITN.GT.0) GO TO 120
      LFAIL = P01AAF(ISAVE,1,SRNAME)
      RETURN
C     FORM SHIFT
  120 IF ((ITS.NE.10) .AND. (ITS.NE.20)) GO TO 140
      SR = DABS(HR(N,N-1)) + DABS(HR(N-1,N-2))
      SI = DABS(HI(N,N-1)) + DABS(HI(N-1,N-2))
      GO TO 180
  140 SR = HR(N,N)
      SI = HI(N,N)
      N1 = N - 1
      XR = HR(N1,N)*HR(N,N1) - HI(N1,N)*HI(N,N1)
      XI = HR(N1,N)*HI(N,N1) + HI(N1,N)*HR(N,N1)
      IF ((XR.EQ.0.0D0)	.AND. (XI.EQ.0.0D0)) GO	TO 180
      YR = (HR(N-1,N-1)-SR)/2.0D0
      YI = (HI(N-1,N-1)-SI)/2.0D0
      CALL A02AAF(YR**2-YI**2+XR, 2.0D0*YR*YI+XI, ZR, ZI)
      IF ((YR*ZR+YI*ZI).GE.0.0D0) GO TO	160
      ZR = -ZR
      ZI = -ZI
  160 CALL A02ACF(XR, XI, YR+ZR, YI+ZI,	ZXR, ZXI)
      XR = ZXR
      XI = ZXI
      SR = SR -	XR
      SI = SI -	XI
  180 DO 200 I=1,N
	 HR(I,I) = HR(I,I) - SR
	 HI(I,I) = HI(I,I) - SI
  200 CONTINUE
      TR = TR +	SR
      TI = TI +	SI
      ITS = ITS	+ 1
      ITN = ITN	- 1
      J	= L + 1
C     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
      XR = DABS(HR(N-1,N-1)) + DABS(HI(N-1,N-1))
      YR = DABS(HR(N,N-1)) + DABS(HI(N,N-1))
      ZR = DABS(HR(N,N)) + DABS(HI(N,N))
      M	= N
      N1 = N - 1
      IF (J.GT.N1) GO TO 240
      DO 220 MM=J,N1
	 M = M - 1
	 YI = YR
	 YR = DABS(HR(M,M-1)) +	DABS(HI(M,M-1))
	 XI = ZR
	 ZR = XR
	 XR = DABS(HR(M-1,M-1))	+ DABS(HI(M-1,M-1))
	 IF (YR.LE.(ACC*ZR/YI*(ZR+XR+XI))) GO TO 260
  220 CONTINUE
  240 M	= L
C     TRIANGULAR DECOMPOSITION H = L*R
  260 M1 = M + 1
      DO 360 I=M1,N
	 XR = HR(I-1,I-1)
	 XI = HI(I-1,I-1)
	 YR = HR(I,I-1)
	 YI = HI(I,I-1)
	 IF ((DABS(XR)+DABS(XI)).GE.(DABS(YR)+DABS(YI))) GO TO 300
C     INTERCHANGE ROWS OF HR AND HI
	 I1 = I	- 1
	 DO 280	J=I1,N
	    ZR = HR(I-1,J)
	    HR(I-1,J) =	HR(I,J)
	    HR(I,J) = ZR
	    ZI = HI(I-1,J)
	    HI(I-1,J) =	HI(I,J)
	    HI(I,J) = ZI
  280	 CONTINUE
	 CALL A02ACF(XR, XI, YR, YI, ZR, ZI)
	 WR(I) = 1.0D0
	 GO TO 320
  300	 CALL A02ACF(YR, YI, XR, XI, ZR, ZI)
	 WR(I) = -1.0D0
  320	 HR(I,I-1) = ZR
	 HI(I,I-1) = ZI
	 DO 340	J=I,N
	    HR(I,J) = HR(I,J) -	ZR*HR(I-1,J) + ZI*HI(I-1,J)
	    HI(I,J) = HI(I,J) -	ZR*HI(I-1,J) - ZI*HR(I-1,J)
  340	 CONTINUE
  360 CONTINUE
C     COMPOSITION R*L =	H
      M1 = M + 1
      DO 440 J=M1,N
	 XR = HR(J,J-1)
	 XI = HI(J,J-1)
	 HR(J,J-1) = 0.0D0
	 HI(J,J-1) = 0.0D0
C     INTERCHANGE COLUMNS OF HR	AND HI,	IF NECESSARY
	 IF (WR(J).LE.0.0D0) GO	TO 400
	 DO 380	I=L,J
	    ZR = HR(I,J-1)
	    HR(I,J-1) =	HR(I,J)
	    HR(I,J) = ZR
	    ZI = HI(I,J-1)
	    HI(I,J-1) =	HI(I,J)
	    HI(I,J) = ZI
  380	 CONTINUE
  400	 DO 420	I=L,J
	    HR(I,J-1) =	HR(I,J-1) + XR*HR(I,J) - XI*HI(I,J)
	    HI(I,J-1) =	HI(I,J-1) + XR*HI(I,J) + XI*HR(I,J)
  420	 CONTINUE
  440 CONTINUE
      GO TO 40
C     A	ROOT FOUND
  460 WR(N) = HR(N,N) +	TR
      WI(N) = HI(N,N) +	TI
      N	= N - 1
      GO TO 20
  480 RETURN
      END

c
      SUBROUTINE A02ACF(XXR, XXI, YYR, YYI, ZR,	ZI)
C     MARK 2A RELEASE.	NAG COPYRIGHT 1973
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C
C     DIVIDES ONE COMPLEX NUMBER BY A SECOND
C
      DOUBLE PRECISION ONE, YYR, YYI, H, A, ZR,	XXR, XXI, ZI
c   * added (Draine!)
      SAVE ONE
      DATA ONE /1.0D0/
C
      IF (DABS(YYR).LE.DABS(YYI)) GO TO	20
      H	= YYI/YYR
      A	= ONE/(H*YYI+YYR)
      ZR = (XXR+H*XXI)*A
      ZI = (XXI-H*XXR)*A
      RETURN
   20 H	= YYR/YYI
      A	= ONE/(H*YYR+YYI)
      ZR = (H*XXR+XXI)*A
      ZI = (H*XXI-XXR)*A
      RETURN
      END

c
      SUBROUTINE F01AMF(N, K, L, AR, IAR, AI, IAI, INTGER)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C
C     COMHES
C     GIVEN THE	COMPLEX	UNSYMMETRIC MATRIX,A, STORED IN	THE ARRAYS
C     AR(N,N) AND AI(N,N), THIS	SUBROUTINE REDUCES THE
C     SUB-MATRIX OF ORDER L - K	+ 1, WHICH STARTS AT THE ELEMENTS
C     AR(K,K) AND AI(K,K) AND FINISHES AT THE ELEMENTS AR(L,L)
C     AND AI(L,L), TO HESSENBERG FORM, H, BY NON-ORTHOGONAL
C     ELEMENTARY TRANSFORMATIONS. THE MATRIX H IS OVERWRITTEN ON A
C     WITH DETAILS OF THE TRANSFORMATIONS STORED IN THE	REMAINING
C     TRIANGLE UNDER H AND IN ELEMENTS K TO L OF THE ARRAY
C     INTGER(N).
C     1ST AUGUST 1971
C
      INTEGER LA, L, K1, K, M, I, J, M1, N, IAR, IAI, INTGER(N), MM
      DOUBLE PRECISION XR, XI, YR, YI, AR(IAR,N), AI(IAI,N)
      LA = L - 1
      K1 = K + 1
      IF (K1.GT.LA) RETURN
      DO 160 M=K1,N
	 I = M
	 XR = 0.0D0
	 XI = 0.0D0
	 IF (M.GT.L) GO	TO 120
	 DO 20 J=M,L
	    IF ((DABS(AR(J,M-1))+DABS(AI(J,M-1))).LE.(DABS(XR)+DABS(XI))
     *	    ) GO TO 20
	    XR = AR(J,M-1)
	    XI = AI(J,M-1)
	    I =	J
   20	 CONTINUE
	 INTGER(M) = I
	 IF (I.EQ.M) GO	TO 80
C     INTERCHANGE ROWS AND COLUMNS OF ARRAYS AR	AND AI.
	 DO 40 J=K,N
	    YR = AR(I,J)
	    AR(I,J) = AR(M,J)
	    AR(M,J) = YR
	    YI = AI(I,J)
	    AI(I,J) = AI(M,J)
	    AI(M,J) = YI
   40	 CONTINUE
	 DO 60 J=1,L
	    YR = AR(J,I)
	    AR(J,I) = AR(J,M)
	    AR(J,M) = YR
	    YI = AI(J,I)
	    AI(J,I) = AI(J,M)
	    AI(J,M) = YI
   60	 CONTINUE
   80	 IF (XR.EQ.0.0D0 .AND. XI.EQ.0.0D0) GO TO 120
	 IF (M.GE.L) GO	TO 120
	 M1 = M	+ 1
	 DO 100	I=M1,L
	    YR = AR(I,M-1)
	    YI = AI(I,M-1)
	    CALL A02ACF(YR, YI,	XR, XI,	AR(I,M-1), AI(I,M-1))
  100	 CONTINUE
	 CALL F01AMZ(AR(1,M+1),	IAR, AI(1,M+1),	IAI, L,	L-M,AR(M+1,M-1),
     *	  1, AI(M+1,M-1), 1, AR(1,M), AI(1,M))
  120	 CALL F01AMY(AR(K+1,K),	IAR, AI(K+1,K),	IAI, L-K, M-K,AR(K+1,M),
     *	  AI(K+1,M))
  160 CONTINUE
      DO 200 MM=K1,LA
	 M = LA	- MM + K1
	 M1 = M	- 2
	 I = INTGER(M)
	 IF (I.EQ.M) GO	TO 200
	 IF (K.GT.M1) GO TO 200
	 DO 180	J=K,M1
	    YR = AR(I,J)
	    AR(I,J) = AR(M,J)
	    AR(M,J) = YR
	    YI = AI(I,J)
	    AI(I,J) = AI(M,J)
	    AI(M,J) = YI
  180	 CONTINUE
  200 CONTINUE
      RETURN
      END

c
      SUBROUTINE F01AMY(AR, IAR, AI, IAI, M, N,	BR, BI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C
C     COMPUTES	B = L**(-1) * B	 (COMPLEX) WHERE
C     L	IS A LOWER TRIANGULAR MATRIX OF	THE SPECIAL FORM
C     (	L11  0 )
C     (	L21  I )
C     L11 IS LOWER TRIANGULAR OF ORDER N,
C     L21 IS RECTANGULAR (M-N) BY N.
C     A	HOLDS THE SUBDIAGONAL ELEMENTS OF THE FIRST N COLUMNS OF L
C     AND THE DIAGONAL ELEMENTS	OF L ARE TAKEN TO BE 1.0 .
C
C     .. SCALAR	ARGUMENTS ..
      INTEGER IAR, IAI,	M, N
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION AR(IAR,N), AI(IAI,N), BR(M), BI(M)
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE PRECISION XR, XI
      INTEGER I, J, JP1
C     ..
C
      DO 40 J=1,N
	 JP1 = J + 1
	 IF (JP1.GT.M) GO TO 40
	 XR = BR(J)
	 XI = BI(J)
	 DO 20 I=JP1,M
	    BR(I) = (BR(I)-AR(I,J)*XR) + AI(I,J)*XI
	    BI(I) = (BI(I)-AR(I,J)*XI) - AI(I,J)*XR
   20	 CONTINUE
   40 CONTINUE
      RETURN
      END

c
      INTEGER FUNCTION P01AAF(IFAIL, ERROR, SRNAME)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED
C     MARK 4A REVISED, IER-45
C     MARK 4.5 REVISED
C     MARK 7 REVISED (DEC 1978)
C     MARK 11 REVISED (FEB 1984)
C     RETURNS THE VALUE	OF ERROR OR TERMINATES THE PROGRAM.
      INTEGER ERROR, IFAIL, NOUT
C$P 1
c * Draine!
      CHARACTER*8 SRNAME
c      DOUBLE PRECISION SRNAME
      CHARACTER REPLAY
C     TEST IF NO ERROR DETECTED
      IF (ERROR.EQ.0) GO TO 20
C     DETERMINE	OUTPUT UNIT FOR	MESSAGE
      CALL X04AAF (0,NOUT)
C     TEST FOR SOFT FAILURE
      IF (MOD(IFAIL,10).EQ.1) GO TO 10
C     HARD FAILURE
****   !!!!!!
****   STDOUT is chanched on '*'
****   V.Grikurov, 15/09/94      
c      WRITE (*,99999) SRNAME, ERROR
c      READ (*, 11111) REPLAY
      WRITE (NOUT,99999) SRNAME, ERROR
      READ (NOUT, 11111) REPLAY
C     ******************** IMPLEMENTATION NOTE ********************
C     THE FOLLOWING STOP STATEMENT MAY BE REPLACED BY A	CALL TO	AN
C     IMPLEMENTATION-DEPENDENT ROUTINE TO DISPLAY A MESSAGE AND/OR
C     ABORT THE	PROGRAM.
C     *************************************************************
      IF(REPLAY.EQ.'Y' .OR. REPLAY.EQ.'y') GOTO 20
      STOP
C     SOFT FAIL
C     TEST IF ERROR MESSAGES SUPPRESSED
   10 IF (MOD(IFAIL/10,10).EQ.0) GO TO 20
      WRITE (NOUT,99999) SRNAME, ERROR
   20 P01AAF = ERROR
      RETURN
99999 FORMAT(1H , 38HERROR DETECTED BY	NAG LIBRARY ROUTINE , A8,
     * 11H - IFAIL = , I5/1H , 22HCONTINUE ANYWAY (Y/N)? )
11111 FORMAT( A1 )
      END

c
      SUBROUTINE A02AAF(XXR, XXI, YR, YI)
C     MARK 2A RELEASE.	NAG COPYRIGHT 1973
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     COMPUTES THE SQUARE ROOT OF A COMPLEX NUMBER
C
      DOUBLE PRECISION ZERO, HALF, XI, XXI, H, XXR, YR,	YI, A02ABF
c   * added Draine!
      SAVE HALF,ZERO
      DATA ZERO	/0.0D0/, HALF /0.5D0/
C
      XI = XXI
      H	= DSQRT((DABS(XXR)+A02ABF(XXR,XI))*HALF)
      IF (XI.NE.ZERO) XI = XI/(H+H)
      IF (XXR.LT.ZERO) GO TO 20
      YR = H
      YI = XI
      RETURN
   20 IF (XI.LT.ZERO) GO TO 40
      YR = XI
      YI = H
      RETURN
   40 YR = -XI
      YI = -H
      RETURN
      END

c
      DOUBLE PRECISION FUNCTION	X02AAF(X)
C     NAG COPYRIGHT 1975
C     MARK 4.5 RELEASE
      DOUBLE PRECISION X
C     *	EPS *
C
C     IBM DOUBLE PRECISION VERSION
C
C     RETURNS THE VALUE	EPS WHERE EPS IS THE SMALLEST
C     POSITIVE
C     NUMBER SUCH THAT 1.0 + EPS > 1.0
C     THE X PARAMETER IS NOT USED
C     FOR ICL 1900
C     X02AAF = 2.0**(-37.0)
C     FOR IBM 360/370
      X02AAF = 2.0D0**(-52.0D0)
c      DOUBLE PRECISION Z
C     DATA Z/Z'3410000000000000'/
c      DATA Z/Z'0000000000000001'/
c      X02AAF = Z
      RETURN
      END

c
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     IF I = 0,	SETS NERR TO CURRENT ERROR MESSAGE UNIT	NUMBER
C     (STORED IN NERR1).
C     IF I = 1,	CHANGES	CURRENT	ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     *** NOTE ***
C     THIS ROUTINE ASSUMES THAT	THE VALUE OF NERR1 IS SAVED
C     BETWEEN CALLS.  IN SOME IMPLEMENTATIONS IT MAY BE
C     NECESSARY	TO STORE NERR1 IN A LABELLED COMMON
C     BLOCK /AX04AA/ TO	ACHIEVE	THIS.
C
C     .. SCALAR	ARGUMENTS ..
      INTEGER I, NERR
C     ..
C     .. LOCAL SCALARS ..
      INTEGER NERR1
C     ..
c * Draine!
      SAVE NERR1
      DATA NERR1 /6/
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1	= NERR
      RETURN
      END

c
      DOUBLE PRECISION FUNCTION	A02ABF(XXR, XXI)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C
C     RETURNS THE ABSOLUTE VALUE OF A COMPLEX NUMBER VIA ROUTINE
C     NAME
C
      DOUBLE PRECISION ZERO, ONE, XR, XXR, XI, XXI, H
c   * added (Draine!)
      SAVE ONE,ZERO
      DATA ZERO	/0.0D0/, ONE /1.0D0/
C
      XR = DABS(XXR)
      XI = DABS(XXI)
      IF (XI.LE.XR) GO TO 20
      H	= XR
      XR = XI
      XI = H
   20 IF (XI.NE.ZERO) GO TO 40
      A02ABF = XR
      RETURN
   40 H	= XR*DSQRT(ONE+(XI/XR)**2)
      A02ABF = H
      RETURN
      END
