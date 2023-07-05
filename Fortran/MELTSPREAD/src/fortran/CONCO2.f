C...................................................................... 
C.....SUBROUTINE CONCO2 EVALUATES THE PROPERTIES FOR CO2 GAS. 
C...................................................................... 
      SUBROUTINE CONCO2(T,PDRYWL,CP,VISC,FK,DEN)    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA RGASAR/0.188919534D07/
      DATA CPAR/1.326500007D03/ 
      TEV=T 
      IF(T.LT.300.) TEV=300.    
      IF(T.GT.2000.) TEV=2000.  
      VISC=3.9573D-8*(TEV-300.D0)+1.496D-5
      FK=8.8426D-5*(TEV-300.D0)+1.657D-2 
      CP=CPAR
      DEN=(1.D10*PDRYWL)/(T*RGASAR)     
      RETURN
      END