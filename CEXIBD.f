      COMPLEX*16 FUNCTION CEXIBD(HW,Z1,Z2,DEPB,RHOB,GRHOB,
     1                          ALPHA1,ALPHA2,BETA1,BETA2)
C
C     THIS FUNCTION EVALUATES THE INTEGRAL OF THE PRODUCT OF
C     THE DERIVATIVES OF TWO BASIC DEPTH FUNCTIONS OVER THE 
C     INTERVAL Z1 TO Z2 IN THE BOTTOM, WEIGHTED BY AN EXPONENTIAL
C     DENSITY PROFILE.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMPLEX*16 ALPHA1,ALPHA2,BETA1,BETA2
      COMPLEX*16 B1,B2,CFAC,CGRHOB,CI
      COMPLEX*16 U,V,W,X
C
      COMMON /BLKEVN/ HB,CW,CB,FKW,FKB,ROHW,ROHB,ATEN
C
      CI=DCMPLX(0.0D0,1.0D0)
      CGRHOB=DCMPLX(GRHOB,0.0D0)
      THBMHW=2.0D0*HB-HW
C
C     THE FOLLOWING INCLUDES A FACTOR OF 1/FKB**2 TO MAKE 
C     THE INTEGRAL NONDIMENSIONAL
C
      CFAC=-FKB**2*(CDSIN(ALPHA1*HW)*CDSIN(ALPHA2*HW)*
     1     BETA1*BETA2)/(RHOB*ALPHA1*ALPHA2)
C
      B1=2.0D0*BETA1*(HB-HW)
      IF(DIMAG(B1) .LT. 88.0D0) THEN
      CFAC=CFAC/(1.0D0-CDEXP(CI*B1))
      END IF
C
      B2=2.0D0*BETA2*(HB-HW)
      IF(DIMAG(B2) .LT. 88.0D0) THEN
      CFAC=CFAC/(1.0D0-CDEXP(CI*B2))
      END IF
C
      CEXIBD=CMPLX(0.0D0,0.0D0)
C
      W=CI*(BETA1+BETA2)-CGRHOB
      IF(CDABS(W) .GT. 1.0D-12) THEN
      U=CDEXP(CI*(BETA1+BETA2)*(Z2-HW)-CGRHOB*(Z2-DEPB))
      V=CDEXP(CI*(BETA1+BETA2)*(Z1-HW)-CGRHOB*(Z1-DEPB))
      CEXIBD=CEXIBD+(U-V)/W
      ELSE 
      X=CDEXP(CI*(BETA1+BETA2)*(DEPB-HW))
      CEXIBD=CEXIBD+X*(Z2-Z1)
      END IF
C
      W=CI*(BETA1-BETA2)-CGRHOB
      IF(CDABS(W) .GT. 1.0D-12) THEN
      U=CDEXP(CI*(BETA1*(Z2-HW)+BETA2*(THBMHW-Z2))-CGRHOB*(Z2-DEPB))
      V=CDEXP(CI*(BETA1*(Z1-HW)+BETA2*(THBMHW-Z1))-CGRHOB*(Z1-DEPB))
      CEXIBD=CEXIBD+(U-V)/W
      ELSE 
      X=CDEXP(CI*(BETA1*(DEPB-HW)+BETA2*(THBMHW-DEPB)))
      CEXIBD=CEXIBD+X*(Z2-Z1)
      END IF
C
      W=CI*(BETA2-BETA1)-CGRHOB
      IF(CDABS(W) .GT. 1.0D-12) THEN
      U=CDEXP(CI*(BETA2*(Z2-HW)+BETA1*(THBMHW-Z2))-CGRHOB*(Z2-DEPB))
      V=CDEXP(CI*(BETA2*(Z1-HW)+BETA1*(THBMHW-Z1))-CGRHOB*(Z1-DEPB))
      CEXIBD=CEXIBD+(U-V)/W
      ELSE 
      X=CDEXP(CI*(BETA2*(DEPB-HW)+BETA1*(THBMHW-DEPB)))
      CEXIBD=CEXIBD+X*(Z2-Z1)
      END IF
C
      W=-CI*(BETA1+BETA2)-CGRHOB
      IF(CDABS(W) .GT. 1.0D-12) THEN
      U=CDEXP(CI*(BETA1+BETA2)*(THBMHW-Z2)-CGRHOB*(Z2-DEPB))
      V=CDEXP(CI*(BETA1+BETA2)*(THBMHW-Z1)-CGRHOB*(Z1-DEPB))
      CEXIBD=CEXIBD+(U-V)/W
      ELSE 
      X=CDEXP(CI*(BETA1+BETA2)*(THBMHW-DEPB))
      CEXIBD=CEXIBD+X*(Z2-Z1)
      END IF
C
      CEXIBD=CFAC*CEXIBD
      RETURN
      END
