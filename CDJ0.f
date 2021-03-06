      COMPLEX*16 FUNCTION CDJ0(Z)
C
C     THIS FUNCTION COMPUTES J0 FOR A SMALL COMPLEX ARGUMENT USING
C     EIGHT TERMS IN THE POWER SERIES EXPANSION ABOUT ZERO 
C     (SEE ABRAMOWITZ & STEGUM PG. 360, EQ. 9.1.10). IT SHOULD BE
C     USEFUL WHEN CDABS(Z) .LT. 2.0D0. THIS IS THE COMPLEX DOUBLE
C     PRECISION VERSION.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,W,ONE,SUM,FACT
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      W=(-0.25D0)*(Z**2)
      FACT=ONE
      SUM=ONE
C
      DO 100 I=1,7
      FACT=FACT*DBLE(FLOAT(I))
      SUM=SUM+(W**I)/(FACT**2)
  100 CONTINUE
C
      CDJ0=SUM
C
      RETURN
      END
