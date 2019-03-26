      COMPLEX*16 FUNCTION CJDJ(EGV,RNG1,RNG2)
C
C     SAME AS CHDH BUT FOR HANKEL FUNCTIONS OF ORDER ZERO TYPE TWO.
C     THE NAME IS CARRIED OVER FROM A VERSION WHICH USED BESSEL
C     FUNCTIONS.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
      COMMON /FLAGS/ IGEOM,NEWSBC
      COMPLEX*16 CDSRT,FKH,ONE,CI,EGV
C
      ONE=DCMPLX(1.0,0.0)
      CI=DCMPLX(0.0,1.0)
C
      FKH=FKB*CDSRT(ONE-EGV)
      IF(IGEOM .EQ. 1) THEN
      CJDJ=CDEXP(CI*FKH*(RNG1-RNG2))
      ELSE
      CJDJ=SQRT(RNG1/RNG2)*CDEXP(CI*FKH*(RNG1-RNG2))
      END IF
C
      RETURN
      END
