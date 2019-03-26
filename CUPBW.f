      COMPLEX*16 FUNCTION CUPBW(H1,H2,ALPHA1,ALPHA2,BETA1,BETA2,
     1                           DEPB,RHOB,GRHOB,LHB,LH1B,MX)
C
C     THIS FUNCTION COMPUTES THE INTEGRAL OF TWO BASIC DEPTH
C     FUNCTIONS OVER THE BOTTOM FROM DEPB(1)=H2 TO HB, WEIGHTED
C     BY THE BOTTOM DENSITY PROFILE RHOB. THIS INTEVAL STRADDLES
C     THE DEPTH H1, WHICH IS THE WATER BOTTOM INTERFACE FOR THE 
C     FIRST DEPTH FUNCITON. THE LOGIC IS DESIGNED TO ACCOMODATE
C     THIS CHANGE AT H1, WHICH IS IN THE INTERVAL BETWEEN DEPB(LH1B) 
C     AND DEPW(LH1B+1). NOTE THAT H1 IS .GT. H2.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DEPB(MX),RHOB(MX),GRHOB(MX)
      COMPLEX*16  CEXIWB,CEXIB
      COMPLEX*16  ALPHA1,ALPHA2,BETA1,BETA2
C
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
C
      LH1BM1=LH1B-1
      LH1BP1=LH1B+1
      LHBM1=LHB-1
C
      CUPBW=DCMPLX(0.0,0.0)
C
      IF(LH1BM1 .LT. 1 ) GO TO 600
      DO 500 L=1,LH1BM1
      Z1=DEPB(L)
      Z2=DEPB(L+1)
      CUPBW=CUPBW+CEXIWB(H2,Z1,Z2,DEPB(L),RHOB(L),GRHOB(L),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
  500 CONTINUE
  600 CONTINUE
C
      Z1=DEPB(LH1B)
      CUPBW=CUPBW+CEXIWB(H2,Z1,H1,DEPB(LH1B),RHOB(LH1B),GRHOB(LH1B),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
C
      IF(LH1B .EQ. LHB) THEN
C
C     H1 IS IN THE LAYER WITH HB
C
      CUPBW=CUPBW+CEXIB(H1,H2,H1,HB,DEPB(LHB),RHOB(LHB),GRHOB(LHB),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
C
      ELSE
C      
      Z2=DEPB(LH1BP1)
      CUPBW=CUPBW+CEXIB(H1,H2,H1,Z2,DEPB(LH1B),RHOB(LH1B),GRHOB(LH1B),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
C
      IF(LHBM1 .LT. LH1BP1) GO TO 800
      DO 700 L=LH1BP1,LHBM1
      Z1=DEPB(L)
      Z2=DEPB(L+1)
      CUPBW=CUPBW+CEXIB(H1,H2,Z1,Z2,DEPB(L),RHOB(L),GRHOB(L),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
  700 CONTINUE
  800 CONTINUE
C
      Z1=DEPB(LHB)
      CUPBW=CUPBW+CEXIB(H1,H2,Z1,HB,DEPB(LHB),RHOB(LHB),GRHOB(LHB),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
C
      END IF
C
      RETURN
      END