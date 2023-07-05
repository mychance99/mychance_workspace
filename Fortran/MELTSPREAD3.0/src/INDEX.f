C...................................................................... 
C.....SUBROUTINE INDEX CONVERTS THE BWRSAR INDICIAL SCHEME TO THAT      
C.....COMPATABLE WITH THE PROPERTY ROUTINES.
C...................................................................... 
      SUBROUTINE INDEX(I,IJIM)  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      GO TO(1,2,2,2,3,2,2,4,5,5,5,6,6,7,8,10),I
    1 CONTINUE
      IJIM=20
      GO TO 10
    2 CONTINUE
      IJIM=I+15
      GO TO 10
    3 CONTINUE
      IJIM=23
      GO TO 10
    4 CONTINUE
      IJIM=27
      GO TO 10
    5 CONTINUE
      IJIM=I+5
      GO TO 10
    6 CONTINUE
      IJIM=I+12
      GO TO 10
    7 CONTINUE
      IJIM=28
      GO TO 10
    8 CONTINUE
      IJIM=26
   10 CONTINUE
      RETURN
      END