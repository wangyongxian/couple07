      SUBROUTINE BOUNDS(DEP,CGWNS,CWNS,NP,MX,H,FMIN,FAVG,FMAX,LH)
C
C     THIS SUBROUTINE COMPUTES THE MINIMUM, MAXIMUM AND THE AVERAGE
C     THE OF REAL PART OF THE PIECEWISE LINEAR WAVE NUMBER SQUARED
C     ON THE INTERVAL 0 TO H. IT ALSO RETURNS THE INDEX LH OF THE
C     INTERVAL IN THE SOUND VELOCITY PROFILE WHICH CONTAINS H. IF
C     H .GT. DEP(NP) THEN LH=NP AND THE REAL PART OF THE WAVE NUMBER
C     SQUARED OF THE LINEARLY EXTRAPOLATED TO H TO OBTAIN THE BOUNDS
C     AND THE AVERAGE.      
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION DEP(MX)
C      REAL*8 DREAL
      COMPLEX*16 CGWNS(MX),CWNS(MX)
C
      FMIN=1.0D38
      FAVG=0.0D0
      FMAX=0.0D0
C
      DO 100 LH=1,NP-1
      IF(DEP(LH+1) .GE. H) GO TO 200
  100 CONTINUE
C
C     IF H .GT. DEP(NP), EXTRAPOLATION WILL BE DONE
C
      LH=NP
  200 CONTINUE
C
      IF(LH .GT. 1)   THEN
C
      DO 300 I=1,LH-1
      DD=DEP(I+1)-DEP(I)
      DDSD2=(DEP(I+1)**2-DEP(I)**2)/2.0
      FAVG=FAVG+DREAL(CGWNS(I))*DDSD2+DREAL(CWNS(I))*DD
C
      FWNS=DREAL(CWNS(I)+CGWNS(I)*DEP(I))
      FMIN=DMIN1(FMIN,FWNS)
      FMAX=DMAX1(FMAX,FWNS)
C
  300 CONTINUE
C
      END IF
C
      FWNS=DREAL(CWNS(LH)+CGWNS(LH)*DEP(LH))
      FMIN=DMIN1(FMIN,FWNS)
      FMAX=DMAX1(FMAX,FWNS)
C
      DD=H-DEP(LH)
      DDSD2=(H**2-DEP(LH)**2)/2.0
      FAVG=FAVG+DREAL(CGWNS(LH))*DDSD2+DREAL(CWNS(LH))*DD
C      
      FAVG=FAVG/(H-dep(1))
C
      FWNS=DREAL(CWNS(LH)+CGWNS(LH)*H)
      FMIN=DMIN1(FMIN,FWNS)
      FMAX=DMAX1(FMAX,FWNS)
C
      RETURN
      END
