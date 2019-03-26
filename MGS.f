      SUBROUTINE MGS(T,Q,U,N,NX)
C
C     THIS SUBROUTINE DOES A MODIFIED GRAM-SCHMIDT DECOMPOSITION
C     T=Q*U WHERE Q IS UNITARY AND U IS UPPER TRIANGULAR. THE 
C     INPUT MATRIX T IS USED AS WORK SPACE AND ITS ORIGINAL VALUES
C     ARE NOT PRESERVED. 
C   
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 T(NX,NX),Q(NX,NX),U(NX,NX),CPROD
      REAL*8 AMP
C
      DO 50 I=1,N
      DO 50 II=1,N
      U(II,I)=DCMPLX(0.0D0,0.0D0)
   50 CONTINUE
C
      U(1,1)=DCMPLX(AMP(T(1,1),N),0.0D0)
      DO 100 II=1,N
      Q(II,1)=T(II,1)/U(1,1)
  100 CONTINUE
C
      NM1=N-1
      DO 500 I=1,NM1
      IP1=I+1
C
      DO 400 J=IP1,N
      U(I,J)=CPROD(T(1,J),Q(1,I),N)
C 
      DO 300 II=1,N
      T(II,J)=T(II,J)-U(I,J)*Q(II,I)
  300 CONTINUE
C
  400 CONTINUE
C
      U(IP1,IP1)=DCMPLX(AMP(T(1,IP1),N),0.0D0)
      DO 450 II=1,N
      Q(II,IP1)=T(II,IP1)/U(IP1,IP1)
  450 CONTINUE
C
  500 CONTINUE
C
      RETURN
      END
C
      COMPLEX*16 FUNCTION CPROD(V1,V2,N)
C  
C     THIS FUNCTION COMPUTES THE COMPLEX SCALAR PRODUCT OF THE 
C     VECTORS V1 AND V2 (V1 DOT CONJG(V2)).
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 V1(N),V2(N)
C
      CPROD=DCMPLX(0.D0,0.D0)
      DO 100 I=1,N
      CPROD=CPROD+V1(I)*DCONJG(V2(I))
  100 CONTINUE
C
      RETURN
      END
C
      REAL*8 FUNCTION AMP(V1,N)
C  
C     THIS FUNCTION COMPUTES THE AMPLITUDE OF THE COMPLEX 
C     VECTORS V1. 
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 V1(N)
C
      AMP=0.D0
      DO 100 I=1,N
      AMP=AMP+CDABS(V1(I))**2
  100 CONTINUE
C
      AMP=DSQRT(AMP)
C
      RETURN
      END