rem Limit input file name to 7 characters to allow space for S, O and I
COPY %1.DAT COUPLE.DAT
ERASE COUPLE.PRT
ERASE COUPLE.LOG
COUPLE.EXE
COPY COUPLE.LOG %1.LOG
COPY COUPLE.PRT %1.PRT
COPY COUPLE.TL %1.TL
COPY COUPLES.TL %1S.TL
COPY COUPLEO.TL %1O.TL
COPY COUPLEI.TL %1I.TL
REM
COPY COUPLE.CPR %1.CPR
COPY COUPLES.CPR %1S.CPR
COPY COUPLEO.CPR %1O.CPR
COPY COUPLEI.CPR %1I.CPR
REM
REM ERASE COUPLE.DAT
REM ERASE COUPLE.PRT
REM ERASE COUPLE.LOG
ERASE COUPLE.TL
ERASE COUPLES.TL
ERASE COUPLEO.TL
ERASE COUPLEI.TL
ERASE COUPLE.CPR
ERASE COUPLES.CPR
ERASE COUPLEO.CPR
ERASE COUPLEI.CPR
REM
