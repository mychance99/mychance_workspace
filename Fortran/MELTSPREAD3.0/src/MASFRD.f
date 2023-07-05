C...................................................................... 
C.....SUBROUTINE MASFRD CALCULATES THE AVERAGE DENSITY GIVEN THE
C.....HEIGHT CONSTITUENT 3-D ARRAY, H.    
C...................................................................... 
      SUBROUTINE MASFRD(ROMLIQ,IX,IY,RMASS,ROC,ROCAVE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/BASPR/ H(999,999,16),HBMOLD(999,999,16)    
      DIMENSION ROMLIQ(28)      
      ROMIN=ROC*RMASS 
      SUM=0.D0
      HT=0.D0
      DO 100 I=1,16   
      HCOMP=DMAX1(H(IX,IY,I),0.D0)
      CALL INDEX(I,IJIM)
      HT=HT+HCOMP     
      REV=ROMLIQ(IJIM)
      IF(I.EQ.16) REV=ROC*RMASS 
      SUM=SUM+HCOMP*REV
  100 CONTINUE
      IF(HT.GT.0.D0) GO TO 500  
      ROCAVE=ROMIN    
      GO TO 510
  500 CONTINUE
      ROCAVE=SUM/HT   
  510 CONTINUE
      RETURN
      END