      COMPLEX*16 FUNCTION HK02(Z)
C
C     THIS FUNCTION COMPUTES THE FIRST TERM IN THE ASYMPTOTIC EXPANSION
C     OF THE HANKEL FUNCTION OF ORDER ZERO OF THE SECOND KIND.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PI
C
      COMPLEX*16 CI,Z,U,V,CDSRT
C
      PI=2.0D0*DACOS(0.0D0)
      CI=DCMPLX(0.0D0,1.0D0)
C
      U=CDSRT(2.0/(PI*Z))
      IF(DIMAG(Z) .GT. 675.0D0) GO TO 50
      V=CDEXP(-CI*Z)*CDEXP(-CI*PI/4.0)
      GO TO 100
   50 CONTINUE
      STOP 'OVERFLOW IN HK02 - MAKE FIRST REGION SMALLER'
  100 CONTINUE
      HK02=U*V
C
      RETURN
      END
