      SUBROUTINE AMAT(HW,M,EGVL,FNORM,A,MX,
     1                DEPW,CGRADW,CINTW,RHOW,GRHOW,NPW,LHW,
     2                DEPB,CGRADB,CINTB,RHOB,GRHOB,NPB,LHB)
C
C     THE FOLLOWING SUBROUTINE COMPUTES THE COMPLEX SYMETRIC MATRIX
C     A NEEDED BY GALRKN. THIS MATRIX IS USED TO FORM PART OF THE
C     GENERALIZED MATRIX EIGENVALUE PROBLEM THAT IS SOLVED TO FIND
C     THE EIGENVALUES AND EIGENFUNCTIONS IN TERMS OF THE BASIC 
C     EIGENVALUES AND EIGENFUNCTIONS. Removed real eigenvalues
c     from the diagonal and corrected error in 4/29/97 version of
c     amat.for. The 4/29/97 version is the one in couple97 on oalib.
c     (R. B. Evans, N. Stonington, CT, 9/25/07)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION DEPW(MX),RHOW(MX),GRHOW(MX)
      DIMENSION DEPB(MX),RHOB(MX),GRHOB(MX)
      COMPLEX*16 CGRADW(MX),CINTW(MX)
      COMPLEX*16 CGRADB(MX),CINTB(MX)
      COMPLEX*16 EGVL(MX),FNORM(MX)
      COMPLEX*16 ONE
      COMPLEX*16 CDSRT,CEXIBD,CEXIB,ZCEXIB
      COMPLEX*16 CEXIWD,CEXIW,ZCEXIW
      COMPLEX*16 EGV1,EGV2,ALPHA1,ALPHA2,BETA1,BETA2
      COMPLEX*16 A(MX,MX),U
C
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
C
      ONE=DCMPLX(1.0,0.0)
      ZERO=0.0D0
      LHWM1=LHW-1
      LHBM1=LHB-1
C
C     STORE THE INTEGRALS IN THE WATER AND SEDIMENT IN THE MATRIX A
C
      DO 700 I=1,M
C
      EGV1=EGVL(I)
      ALPHA1=CDSRT(FKW**2+FKB**2*(EGV1-ONE))
      BETA1=CDSRT(FKB**2*EGV1)
C
      DO 500 II=I,M
C
      EGV2=EGVL(II)
      ALPHA2=CDSRT(FKW**2+FKB**2*(EGV2-ONE))
      BETA2=CDSRT(FKB**2*EGV2)
C
      U=CDSRT(FNORM(I))*CDSRT(FNORM(II))
C
C     INITIALIZE WITH ZERO
C
      A(II,I)=DCMPLX(0.0D0,0.0D0)
C
C     SUM OVER LAYERS IN THE WATER SOUND SPEED-DENSITY PROFILE
C
      IF(LHW .EQ. 1) GO TO 200
      DO 150 L=1,LHWM1
      Z1=DEPW(L)
      Z2=DEPW(L+1)
      A(II,I)=A(II,I)+CEXIWD(Z1,Z2,DEPW(L),RHOW(L),GRHOW(L),
     1                            ALPHA1,ALPHA2)
      A(II,I)=A(II,I)+(ONE-CINTW(L))*CEXIW(Z1,Z2,DEPW(L),RHOW(L),
     1                   GRHOW(L),ALPHA1,ALPHA2)
     2               -CGRADW(L)*ZCEXIW(Z1,Z2,DEPW(L),RHOW(L),
     3                   GRHOW(L),ALPHA1,ALPHA2)
  150 CONTINUE
  200 CONTINUE
C
      Z1=DEPW(LHW)
      Z2=HW
      A(II,I)=A(II,I)+CEXIWD(Z1,Z2,DEPW(LHW),RHOW(LHW),
     1                            GRHOW(LHW),ALPHA1,ALPHA2)
      A(II,I)=A(II,I)+(ONE-CINTW(LHW))*CEXIW(Z1,Z2,DEPW(LHW),
     1                   RHOW(LHW),GRHOW(LHW),ALPHA1,ALPHA2)
     2               -CGRADW(LHW)*ZCEXIW(Z1,Z2,DEPW(LHW),
     3                   RHOW(LHW),GRHOW(LHW),ALPHA1,ALPHA2)
C
C     SUM OVER LAYERS IN BOTTOM SOUND SPEED - DENSITY PROFILE
C
      IF(LHB .EQ. 1) GO TO 300
      DO 250 L=1,LHBM1
      Z1=DEPB(L)
      Z2=DEPB(L+1)
      A(II,I)=A(II,I)+CEXIBD(HW,Z1,Z2,DEPB(L),RHOB(L),GRHOB(L),
     1                            ALPHA1,ALPHA2,BETA1,BETA2)
      A(II,I)=A(II,I)+(ONE-CINTB(L))*CEXIB(HW,HW,Z1,Z2,DEPB(L),
     1                  RHOB(L),GRHOB(L),ALPHA1,ALPHA2,BETA1,BETA2)
     2               -CGRADB(L)*ZCEXIB(HW,Z1,Z2,DEPB(L),
     3                  RHOB(L),GRHOB(L),ALPHA1,ALPHA2,BETA1,BETA2)
  250 CONTINUE
  300 CONTINUE
C
      Z1=DEPB(LHB)
      Z2=HB      
      A(II,I)=A(II,I)+CEXIBD(HW,Z1,Z2,DEPB(LHB),RHOB(LHB),GRHOB(LHB),
     1                            ALPHA1,ALPHA2,BETA1,BETA2)
      A(II,I)=A(II,I)+(ONE-CINTB(LHB))*CEXIB(HW,HW,Z1,Z2,DEPB(LHB),
     1                 RHOB(LHB),GRHOB(LHB),ALPHA1,ALPHA2,BETA1,BETA2)
     2               -CGRADB(LHB)*ZCEXIB(HW,Z1,Z2,DEPB(LHB),
     3                 RHOB(LHB),GRHOB(LHB),ALPHA1,ALPHA2,BETA1,BETA2)
C
C     APPLY NORMALIZATION FACTORS OF BASIC DEPTH FUNCTIONS
C
      A(II,I)=A(II,I)/U
C
C     DEFINE THE REST OF A BY SYMMETRY
C
      IF(II .NE. I) A(I,II)=A(II,I)
C
  500 CONTINUE
  700 CONTINUE
C
      RETURN
      END