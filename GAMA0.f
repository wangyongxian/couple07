      COMPLEX*16 FUNCTION GAMA0(H,ALPHA,BETA)
C
C     THIS FUNCTION COMPUTES THE NORMALIZATION FACTOR FOR THE EIGEN-
C     FUNCTION CORRESPONDING TO THE VERTICAL WAVE NUMBERS ALPHA AND
C     BETA.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 ALPHA,BETA,U,V,W,Z,B
      COMPLEX*16 CI,ONE
      COMPLEX*16 CDSRT,CDTAN
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
C
      ONE=DCMPLX(1.0D0,0.0D0)
      CI=DCMPLX(0.0D0,1.0D0)
C
      B=BETA*(H-HB)
      Z=FKB**2/ALPHA
      IF(DIMAG(B) .LT. -44.0D0) GO TO 100
      V=(B/(CDSIN(B)**2)-1./CDTAN(B))/(2.0*BETA*ROHB)
      GO TO 200
  100 CONTINUE
      V=(-1./CDTAN(B))/(2.0*BETA*ROHB)
  200 CONTINUE
      U=(ALPHA*H-CDCOS(ALPHA*H)*CDSIN(ALPHA*H))/(2.0*ALPHA*ROHW)
      W= CDSIN(ALPHA*H)
      GAMA0=Z**2*(U-W**2*V)
C
      RETURN
      END
