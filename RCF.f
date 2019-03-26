      REAL*8 FUNCTION RCF(H,X)
C
C     THIS FUNCTION EVALUATES THE CHARACTERISTIC FUNCTION FOR THE
C     PROBLEM WITH ATEN=0.0 USING DOUBLE PRECISION.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ALPHA,BETA,X,ZERO,SDALP,CALP
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
C
C     IT IS ASSUMED THAT X IS LARGE ENOUGH TO MAKE ALPHA REAL. BETA IS
C     DEFINED SO AS TO ELIMINATE THE NEED FOR COMPLEX TRIG. FUNCTIONS.
C
      ZERO=0.0D0
      ALPHA=DSQRT(DABS(FKW**2+FKB**2*(X-1)))
      BETA=DSQRT(DABS(FKB**2*X))
C
C     THE FOLLOWING LOGIC ACCOMODATES SMALL VALUES OF ALPHA
C
      IF(ALPHA .LT. .1D-28) THEN
      SDALP=H
      ELSE
      SDALP=DSIN(ALPHA*H)/ALPHA
      END IF
      CALP=DCOS(ALPHA*H)
C
C     THE VALUE OF THE CHARACTERISTIC FUNCTION WHEN X=ZERO IS STORED 
C     IN RCF
C
      RCF=ROHW*SDALP-ROHB*CALP*(H-HB)
C
C     WHEN X IS .LT. ZERO THE VALUE OF THE CHARACTERISTIC FUNCTION IS
C     DIVIDED BY COSH(BETA*(H-HB) FOR NUMERICAL PURPOSES. THE ZEROES
C     OF THE RESULTING FUNCTION ARE STILL THE SAME.
C
      IF(X .LT. ZERO)
     1 RCF=ROHW*SDALP
     2    -ROHB*CALP*DTANH(BETA*(H-HB))/BETA
C
C
      IF(X .GT. ZERO)
     1 RCF=ROHW*SDALP*DCOS(BETA*(H-HB))
     2    -ROHB*CALP*DSIN(BETA*(H-HB))/BETA
C
      RETURN
      END
