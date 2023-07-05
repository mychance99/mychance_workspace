C...................................................................... 
C.....GIVEN THE SPECIFIC ENTHALPY, SUBROUTINE TESS FINDS THE  
C.....TEMPERATURE FOR CARBON STEEL.
C...................................................................... 
      SUBROUTINE TESS(E,T,DTDE) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DATSH/TS,TL,ES,EL,DFUS,CL    
      DATA TREF/298.15D0/
      IF(E.GT.ES) GO TO 10      
      T=TREF+(E*(TS-TREF))/ES   
      DTDE=(TS-TREF)/ES
      GO TO 200
   10 CONTINUE
      IF(E.GT.EL) GO TO 20      
      T=TS+((E-ES)*(TL-TS))/DFUS
      DTDE=(TL-TS)/DFUS
      GO TO 200
   20 CONTINUE
      T=TL+(E-EL)/CL  
      DTDE=1.D0/CL    
  200 CONTINUE
      RETURN
      END