      SUBROUTINE FORESUB(B,Y,P,M,MX)
C
C     THIS SUBROUTINE SOLVES THE MATRIX EQUATION B*Y=P BY 
C     FORWARD SUBSTITUTION WHERE B IS A REAL LOWER TRIANGULAR
C     MATRIX AND P IS COMPLEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 B(MX,MX)
      COMPLEX*16 Y(MX,MX),P(MX,MX),SUM
C
      DO 1000 J=1,M
C
C     SOLVE FOR THE JTH COLUMN IN Y USING THE JTH COLUMN IN P
C
       Y(1,J)=P(1,J)/DCMPLX(B(1,1),0.0D0)
       DO 800 I=2,M
        IM1=I-1
        SUM=DCMPLX(0.0,0.0)
        DO 600 II=1,IM1
         SUM=SUM+DCMPLX(B(I,II),0.0D0)*Y(II,J)
 600   CONTINUE
        Y(I,J)=(P(I,J)-SUM)/DCMPLX(B(I,I),0.0D0)
 800   CONTINUE
C
 1000 CONTINUE
C
      RETURN
      END
