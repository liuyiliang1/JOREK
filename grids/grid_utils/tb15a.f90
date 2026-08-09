SUBROUTINE TB15A(N,X,F,D,W,LP)
implicit none
!------------------------------------------------------------------
! HSL routine for cubic spline with periodic boundary conditions
! first point must be the same as last : f(1)=f(n)
!    N : number of points
!    X : coordinate (input)
!    F : the function values to be splined (input)
!    D : the derivatives at the points (output)
!    W : workspace (dimension 3N)
!   LP : unit number for output
!------------------------------------------------------------------
REAL*8     :: ZERO,ONE,TWO,THREE
PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)
INTEGER    :: LP,N
REAL*8     :: D(N),F(N),W(*),X(N)
REAL*8     :: A3N1,F1,F2,H1,H2,P
INTEGER    :: I,J,K,N2

IF (N.LT.4) THEN
  WRITE (LP,'(A39)')  'RETURN FROM TB15AD BECAUSE N TOO SMALL'
  W(1) = ONE
  RETURN
END IF
DO I = 2,N
  IF (X(I).LE.X(I-1)) THEN
    WRITE (LP,'(A29,I3,A13)') ' RETURN FROM TB15AD BECAUSE  ',I,' OUT OF ORDER'
    W(1) = TWO
    RETURN
  END IF
ENDDO
IF (F(1).NE.F(N)) THEN
  WRITE (LP,'(A40)')  'RETURN FROM TB15AD BECAUSE F(1).NE.F(N)'
  W(1) = THREE
  RETURN
END IF
DO I = 2,N
  H1 = ONE/ (X(I)-X(I-1))
  F1 = F(I-1)
  IF (I.EQ.N) THEN
    H2 = ONE/ (X(2)-X(1))
    F2 = F(2)
  ELSE
    H2 = ONE/ (X(I+1)-X(I))
    F2 = F(I+1)
  END IF
  W(3*I-2) = H1
  W(3*I-1) = TWO* (H1+H2)
  W(3*I) = H2
  D(I) = 3.0* (F2*H2*H2+F(I)* (H1*H1-H2*H2)-F1*H1*H1)
ENDDO
N2 = N - 2
K = 5
A3N1 = W(3*N-1)
DO I = 2,N2
  P = W(K+2)/W(K)
  W(K+3) = W(K+3) - P*W(K+1)
  D(I+1) = D(I+1) - P*D(I)
  W(K+2) = -P*W(K-1)
  P = W(K-1)/W(K)
  A3N1 = -P*W(K-1) + A3N1
  D(N) = D(N) - P*D(I)
  K = K + 3
ENDDO
P = (W(K+2)+W(K-1))/W(K)
A3N1 = A3N1 - P* (W(K+1)+W(K-1))
D(N) = (D(N)-P*D(N-1))/A3N1
DO I = 3,N
  J = N + 2 - I
  D(J) = (D(J)-W(3*J)*D(J+1)-W(3*J-2)*D(N))/W(3*J-1)
ENDDO
D(1) = D(N)
W(1) = ZERO
RETURN
END
