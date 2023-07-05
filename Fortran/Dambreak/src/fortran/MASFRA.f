C...................................................................... 
C.....SUBROUTINE MASFRA EVALUATES THE AVERAGE DENSITY GIVEN THE ARRAY   
C.....OF CONSTITUENT HEIGHTS, H.
C...................................................................... 
      SUBROUTINE MASFRA(H,ROMLIQ,JSEND,RMASS,ROC,ROCAVE,NCRAB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION H(16,999),ROMLIQ(28)
      ROMIN=ROC*RMASS 
      IF(NCRAB.GT.0) GO TO 10   
      KLOW=1
      KHIGH=16
      GO TO 30
   10 CONTINUE
      IF(NCRAB.EQ.2) GO TO 20   
      KLOW=8
      KHIGH=16
      GO TO 30
   20 CONTINUE
      KLOW=1
      KHIGH=7
   30 CONTINUE
      SUM=0.D0
      HT=0.D0
      DO 40 I=KLOW,KHIGH
      HCOMP=DMAX1(H(I,JSEND),0.D0)
      CALL INDEX(I,IJIM)
      HT=HT+HCOMP     
      REV=ROMLIQ(IJIM)
      IF(I.EQ.16) REV=ROC*RMASS 
      SUM=SUM+HCOMP*REV
   40 CONTINUE
      IF(HT.GT.0.D0) GO TO 60   
      ROCAVE=ROMIN    
      GO TO 70
   60 CONTINUE
      ROCAVE=SUM/HT   
   70 CONTINUE
      RETURN
      END