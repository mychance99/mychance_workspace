C...................................................................... 
C.....SUBROUTINE TTODAB CALCULATES THE BINARY DIFFUSION COEFFICIENTS    
C.....FOR H2O/H2 AND CO2/CO USING HIRSCHFELDER'S EQUATION.    
C...................................................................... 
      SUBROUTINE TTODAB(TEMP1,TEMP2,D1,D2,P)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      X1=TEMP1/((59.7*809.1)**.5)
      OMEGA1=1.579-0.352*X1+.0567*X1*X1-.00317*X1*X1*X1
      D1=.001858D-4*TEMP1**1.5*(1./18.+1./2.)**.5   
      D1=D1/(P*10.*OMEGA1*(.5*(2.827+2.641))**2)    
      X2=TEMP2/((91.7*195.2)**.5)
      OMEGA2=1.070-.0626*X2+.0039*X2*X2-.87D-4*X2*X2*X2
      D2=.001858D-4*TEMP2**1.5*(1./28.+1./44.)**.5  
      D2=D2/(P*10.*OMEGA2*(.5*(3.690+3.941))**2)    
      RETURN
      END