C...................................................................... 
C.....SUBROUTINE VOIDFR EVALUATES THE MELT VOID FRACTION USING THE      
C.....CORRELATION OF KATAOKA AND ISHII.   
C...................................................................... 
      SUBROUTINE VOIDFR(SM,PM,UM,PG,VG,VOID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      GRAV=9.82D0     
      DP=PM-PG
      Y1=SQRT(SM/(GRAV*DP))     
      XNUF=UM/SQRT(SM*PM*Y1)    
      VGO=((SM*GRAV*DP)/PM**2)**0.25D0    
      GAMG=PG/PM      
      CO=1.2D0-0.2D0*SQRT(GAMG) 
      VGSTAR=VG/VGO   
      IF(VGSTAR.GT.0.D0) GO TO 50
      VOID=0.D0
      GO TO 100
   50 CONTINUE
      IF(VGSTAR.GT.0.5D0) GO TO 10
C.....FLOW IS BUBBLY; FIND VOID FRACTION. 
      D2=SQRT(2.D0)   
      VOID=1.D0/(CO+D2/VGSTAR)  
      GO TO 100
   10 CONTINUE
C.....FLOW IS CHURN TURBULENT; CALCULATE VOID FRACTION.
      Z1=1.D0/GAMG**0.157D0     
      Z2=1.D0/XNUF**0.562D0     
      IF(XNUF.GT.2.25D-3) GO TO 35
      Z3=0.03D0*Z1*Z2/VGSTAR    
      VOID=1.D0/(CO+Z3)
      GO TO 100
   35 CONTINUE
      Z4=CO+0.92D0*Z1/VGSTAR    
      VOID=1.D0/Z4    
  100 CONTINUE
      RETURN
      END   