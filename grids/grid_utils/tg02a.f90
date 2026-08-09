SUBROUTINE TG02A(IX,N,U,S,D,X,V)
implicit none
!------------------------------------------------------------------
! HSL subroutine to calculate splined values
!    N  : number of points
!    IX : negative 0 -> no initial guess for where xi is
!        positive -> gues for index close to value X
!   U   : the coordinates of the spline points
!   S   : the function values of the spline points
!   D   : the derivatives on the spline points
!   X   : the coordinate where the output is wanted
!   V(1-4) : value and derivatives of the spline interpolation
!------------------------------------------------------------------
REAL*8   :: X
INTEGER  :: IX,N
REAL*8   :: D(*),S(*),U(*),V(*)
REAL*8   :: A,B,C,C3,D0,D1,EPS,GAMA,H,HR,HRR,PHI,S0,S1,T,THETA
INTEGER  :: I,IFLG,J,K

EPS = 1.D-33
J = 0
K = 0
IFLG = 0
IF (X.LT.U(1)) GO TO 990
IF (X.GT.U(N)) GO TO 991
IF (IX.LT.0 .OR. IFLG.EQ.0) GO TO 12
IF (X.GT.U(J+1)) GO TO 1
IF (X.GE.U(J)) GO TO 18
GO TO 2

    1 J = J + 1
   11 IF (X.GT.U(J+1)) GO TO 1
      GO TO 7
   12 J = ABS(X-U(1))/ (U(N)-U(1))* (N-1) + 1
      J = MIN(J,N-1)
      IFLG = 1
      IF (X.GE.U(J)) GO TO 11
    2 J = J - 1
      IF (X.LT.U(J)) GO TO 2
    7 K = J
      H = U(J+1) - U(J)
      HR = 1./H
      HRR = (HR+HR)*HR
      S0 = S(J)
      S1 = S(J+1)
      D0 = D(J)
      D1 = D(J+1)
      A = S1 - S0
      B = A - H*D1
      A = A - H*D0
      C = A + B
      C3 = C*3.
   18 THETA = (X-U(J))*HR
      PHI = 1. - THETA
      T = THETA*PHI
      GAMA = THETA*B - PHI*A
      V(1) = THETA*S1 + PHI*S0 + T*GAMA
      V(2) = THETA*D1 + PHI*D0 + T*C3*HR
      V(3) = (C* (PHI-THETA)-GAMA)*HRR
      V(4) = -C3*HRR*HR
      RETURN
  990 IF (X.LE.U(1)-EPS*MAX(DABS(U(1)),DABS(U(N)))) GO TO 99
      J = 1
      GO TO 7
  991 IF (X.GE.U(N)+EPS*MAX(DABS(U(1)),DABS(U(N)))) GO TO 995
      J = N - 1
      GO TO 7
  995 K = N
   99 IFLG = 0
DO I = 1,4
  V(I) = 0.
ENDDO
RETURN
END
