      SUBROUTINE FILE(BSOUT,MX,MXT2)
C
C     THIS SUBROUTINE OPENS THE SCRATCH AND OUTPUT FILES GENERATED BY 
C     COUPLE 
C 
C     ELIMINATED STATUS='SCRATCH' FOR LF90 ON 8/20/96. THESE 
C     FILES ARE CLOSED WITH STATUS='DELETE' IN THE MAIN 
C     PROGRAM COUPLE. BROKE UP ALL 2M X 2M FILES INTO M X M FILES
C     ON 2/1/07
C
C     RECL=16*MX*MX IS NEEDED TO WRITE A MX*MX COMPLEX MATRIX 
C     (2*MX*MX 8 BYTE (64 BIT) WORDS)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER BSOUT
C
C      NRECL=16*MX*MX
C      OPEN(8,STATUS='SCRATCH',RECL=NRECL,FORM='UNFORMATTED')
C
      OPEN(8, FILE='SCRACH8.DAT', STATUS='UNKNOWN',FORM='UNFORMATTED'
     1       ,IOSTAT=I8)
      OPEN(23, FILE='SCRAC23.DAT', STATUS='UNKNOWN',FORM='UNFORMATTED'
     1       ,IOSTAT=I23)
      OPEN(24, FILE='SCRAC24.DAT', STATUS='UNKNOWN',FORM='UNFORMATTED'
     1       ,IOSTAT=I24)
      OPEN(31, FILE='SCRAC31.DAT', STATUS='UNKNOWN',FORM='UNFORMATTED'
     1       ,IOSTAT=I31)
      OPEN(33, FILE='SCRAC33.DAT', STATUS='UNKNOWN',FORM='UNFORMATTED'
     1       ,IOSTAT=I33)
      OPEN(34, FILE='SCRAC34.DAT', STATUS='UNKNOWN',FORM='UNFORMATTED'
     1       ,IOSTAT=I34)
C
c      write(*,*) 'i8= ',i8
c      write(*,*) 'i23= ',i23
c      write(*,*) 'i24= ',i24
c      write(*,*) 'i31= ',i31
c      write(*,*) 'i33= ',i33
c      write(*,*) 'i34= ',i34
C
C      NRECL=16*MXT2*MXT2
C      OPEN(4,STATUS='SCRATCH',RECL=NRECL,FORM='UNFORMATTED')
C      OPEN(4,STATUS='SCRATCH',FORM='UNFORMATTED',IOSTAT=I4)
C
      OPEN(44 ,FILE='SCRAC44.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED'
     1      ,IOSTAT=I44)
      OPEN(43 ,FILE='SCRAC43.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED'
     1      ,IOSTAT=I43)
      OPEN(42 ,FILE='SCRAC42.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED'
     1      ,IOSTAT=I42)
      OPEN(41 ,FILE='SCRAC41.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED'
     1      ,IOSTAT=I41)
C
c      write(*,*) 'i44= ',i44
c      write(*,*) 'i43= ',i43
c      write(*,*) 'i42= ',i42
c      write(*,*) 'i41= ',i41
C
      OPEN(10,FILE='COUPLE.TL',FORM='FORMATTED',IOSTAT=I10)
C
c      write(*,*) 'i10= ',i10
C
      OPEN(7 ,FILE='COUPLE.CPR',FORM='UNFORMATTED',IOSTAT=I7)
C
c      write(*,*) 'i7= ',i7
C
c      write(*,*) 'BSOUT= ',BSOUT
      IF(BSOUT .eq. 0) THEN 
      OPEN(12,FILE='COUPLEO.TL',FORM='FORMATTED',IOSTAT=I12)
      OPEN(13,FILE='COUPLEI.TL',FORM='FORMATTED',IOSTAT=I13)
      OPEN(17,FILE='COUPLES.TL',FORM='FORMATTED',IOSTAT=I17)
C
c      write(*,*) 'i12= ',i12
c      write(*,*) 'i13= ',i13
c      write(*,*) 'i17= ',i17
C
      OPEN(14,FILE='COUPLEO.CPR',FORM='UNFORMATTED',IOSTAT=I14)
      OPEN(15,FILE='COUPLEI.CPR',FORM='UNFORMATTED',IOSTAT=I15)
      OPEN(16,FILE='COUPLES.CPR',FORM='UNFORMATTED',IOSTAT=I16)
C
c      write(*,*) 'i14= ',i14
c      write(*,*) 'i15= ',i15
c      write(*,*) 'i16= ',i16
C
      END IF
C
      RETURN
      END
