      COMPLEX*16 FUNCTION CUPWB(H1,H2,ALPHA1,ALPHA2,BETA1,BETA2,
     1                           DEPW,RHOW,GRHOW,LH1,LH2,MX)
C
C     THIS FUNCTION COMPUTES THE INTEGRAL OF TWO BASIC DEPTH
C     FUNCTIONS OVER THE WATER FROM DEPW(1)=0 TO H1, WEIGHTED
C     BY THE WATER DENSITY PROFILE RHOW. THIS INTEVAL STRADDLES
C     THE DEPTH H2, WHICH IS THE WATER BOTTOM INTERFACE FOR THE 
C     SECOND DEPTH FUNCITON. THE LOGIC IS DESIGNED TO ACCOMODATE
C     THIS CHANGE AT H2, WHICH IS IN THE INTERVAL BETWEEN DEPW(LH2) 
C     AND DEPW(LH2). NOTE THAT H2 IS .LT. H1.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DEPW(MX),RHOW(MX),GRHOW(MX)
      COMPLEX*16  CEXIW,CEXIWB
      COMPLEX*16  ALPHA1,ALPHA2,BETA1,BETA2
C
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
C
      LH2M1=LH2-1
      LH2P1=LH2+1
      LH1M1=LH1-1
C
      CUPWB=DCMPLX(0.0,0.0)
C
      IF(LH2M1 .LT. 1 ) GO TO 200
      DO 100 L=1,LH2M1
      Z1=DEPW(L)
      Z2=DEPW(L+1)
      CUPWB=CUPWB+CEXIW(Z1,Z2,DEPW(L),RHOW(L),GRHOW(L),
     1                  ALPHA1,ALPHA2)
  100 CONTINUE
  200 CONTINUE
C
      Z1=DEPW(LH2)
      CUPWB=CUPWB+CEXIW(Z1,H2,DEPW(LH2),RHOW(LH2),GRHOW(LH2),
     1                  ALPHA1,ALPHA2)
C
      IF(LH2 .EQ. LH1) THEN
C
C     H2 IS IN THE LAYER WITH H1
C
      CUPWB=CUPWB+CEXIWB(H2,H2,H1,DEPW(LH1),RHOW(LH1),GRHOW(LH1),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
C
      ELSE
C      
      Z2=DEPW(LH2P1)
      CUPWB=CUPWB+CEXIWB(H2,H2,Z2,DEPW(LH2),RHOW(LH2),GRHOW(LH2),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
C
      IF(LH1M1 .LT. LH2P1) GO TO 400
      DO 300 L=LH2P1,LH1M1
      Z1=DEPW(L)
      Z2=DEPW(L+1)
      CUPWB=CUPWB+CEXIWB(H2,Z1,Z2,DEPW(L),RHOW(L),GRHOW(L),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
  300 CONTINUE
  400 CONTINUE
C
      Z1=DEPW(LH1)
      CUPWB=CUPWB+CEXIWB(H2,Z1,H1,DEPW(LH1),RHOW(LH1),GRHOW(LH1),
     1          ALPHA1,ALPHA2,BETA1,BETA2)
C
      END IF
C
      RETURN
      END
