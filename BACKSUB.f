      SUBROUTINE BACKSUB(B,Y,P,M,MX)
C
C     THIS SUBROUTINE SOLVES THE MATRIX EQUATION B*Y=P BY 
C     BACKSUBSTITUTION WHERE B IS UPPER TRIANGULAR
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 B(MX,MX),Y(MX,MX),P(MX,MX),SUM
C

      MM1=M-1
      DO 1000 J=1,M
       Y(M,J)=P(M,J)/B(M,M)
       DO 800 I=MM1,1,-1
        IP1=I+1
        SUM=DCMPLX(0.0,0.0)
        DO 600 II=IP1,M
         SUM=SUM+B(I,II)*Y(II,J)
 600   CONTINUE
        Y(I,J)=(P(I,J)-SUM)/B(I,I)
 800   CONTINUE
 1000 CONTINUE
C
      RETURN
      END
