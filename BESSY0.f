      FUNCTION BESSY0(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *    -.2073370639D-5,.2093887211D-6/, Q1,Q2,Q3,Q4,Q5/-.1562499995D-
     *1,
     *    .1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
      DATA R1,R2,R3,R4,R5,R6/-2957821389.D0,7062834065.D0,-512359803.6D0
     *,
     *    10879881.29D0,-86327.92757D0,228.4622733D0/,
     *    S1,S2,S3,S4,S5,S6/40076544269.D0,745249964.8D0,
     *    7189466.438D0,47447.26470D0,226.1030244D0,1.D0/
      IF(X.LT.8.)THEN
        Y=X**2
        BESSY0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))/(S1+Y*(S2+Y
     *      *(S3+Y*(S4+Y*(S5+Y*S6)))))+.636619772*BESSJ0(X)*LOG(X)
      ELSE
        Z=8./X
        Y=Z**2
        XX=X-.785398164
        BESSY0=SQRT(.636619772/X)*(SIN(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*
     *      P5))))+Z*COS(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF
      RETURN
      END
