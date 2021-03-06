      COMPLEX*16 FUNCTION CUPW(HW,ALPHA1,ALPHA2,
     1                         DEPW,RHOW,GRHOW,LHW,MX)
C
C     THIS FUNCTION COMPUTES THE INTEGRAL OF TWO BASIC DEPTH
C     FUNCTIONS OVER THE WATER FROM DEPW(1)=0 TO HW, WEIGHTED
C     BY THE WATER DENSITY PROFILE RHOW.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DEPW(MX),RHOW(MX),GRHOW(MX)
      COMPLEX*16  CEXIW
      COMPLEX*16  ALPHA1,ALPHA2
C
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
C
      LHWM1=LHW-1
C
      CUPW=DCMPLX(0.0,0.0)
C     
      IF(LHWM1 .LT. 1 ) GO TO 320
      DO 300 L=1,LHWM1
      Z1=DEPW(L)
      Z2=DEPW(L+1)
      CUPW=CUPW+CEXIW(Z1,Z2,DEPW(L),RHOW(L),GRHOW(L),
     1          ALPHA1,ALPHA2)
  300 CONTINUE
  320 CONTINUE
      Z1=DEPW(LHW)
      CUPW=CUPW+CEXIW(Z1,HW,DEPW(LHW),RHOW(LHW),GRHOW(LHW),
     1          ALPHA1,ALPHA2)
C
C
      RETURN
      END
