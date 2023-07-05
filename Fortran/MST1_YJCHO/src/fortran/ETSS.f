C...................................................................... 
C.....GIVEN THE TEMPERATURE, SUBROUTINE ETSS FINDS THE SPECIFIC
C.....ENTHALPY FOR CARBON STEEL.
C...................................................................... 
      SUBROUTINE ETSS(E,T,DTDE) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DATSH/TS,TL,ES,EL,DFUS,CL    
      DATA TREF/298.15D0/
      IF(T.GT.TS) GO TO 10      
      E=(ES*(T-TREF))/(TS-TREF) 
      DTDE=(TS-TREF)/ES
      GO TO 200
   10 CONTINUE
      IF(T.GT.TL) GO TO 20      
      E=ES+(DFUS*(T-TS))/(TL-TS)
      DTDE=(TL-TS)/DFUS
      GO TO 200
   20 CONTINUE
      E=EL+CL*(T-TL)  
      DTDE=1.D0/CL    
  200 CONTINUE
      RETURN
      END