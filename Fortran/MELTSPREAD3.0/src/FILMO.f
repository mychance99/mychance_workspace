C...................................................................... 
C.....SUBROUTINE FILMO EVALUATES FILM THICKNESS WHEN NO NOCONDENSABLE   
C.....GAS IS PRESENT IN THE FILM.
C...................................................................... 
      SUBROUTINE FILMO(ALPHA,DEL,ASTAR)   
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      E6427=64.D0/27.D0
      ALP=ALPHA/ASTAR 
      ALPR=ABS(ALP)   
      E13=1.D0/3.D0   
      W=SQRT(.25D0*ALPR*ALPR*ALPR*ALPR+E6427)
      Y1P=(W+.5D0*ALPR*ALPR)**E13
      Y1MARG=W-0.5D0*ALPR*ALPR  
      Y1M=0.D0
      IF(Y1MARG.GT.0.D0) Y1M=Y1MARG**E13  
      Y1=Y1P-Y1M      
      SIGN=ALP/ABS(ALP)
      C1=2.D0*SQRT(Y1*Y1+4.D0)  
      C2=SQRT(C1-Y1)  
      DEL=(.5D0*(C2-SIGN*SQRT(Y1)))*ASTAR**0.25D0   
      RETURN
      END