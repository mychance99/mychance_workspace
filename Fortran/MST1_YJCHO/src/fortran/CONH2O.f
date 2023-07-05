C...................................................................... 
C.....SUBROUTINE CONH2O EVALUATES THE PROPERTIES FOR H2O VAPOR.
C...................................................................... 
      SUBROUTINE CONH2O(TEV,PDRYWL,CP,VISC,FK,DEN)    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      T=DMAX1(TEV,293.D0)
      T=DMIN1(T,443.D0) 
      T2=T*T
      T3=T*T2
      VISC=1.D-06*(-0.36792+1.62815D-02*T+7.70341D-05*T2-8.39984D-08*T3)
      DEN=-56.33799+0.53218*T-1.68066D-03*T2+1.77744D-06*T3      
      CP=-817.0310+27.59825*T-9.6764D-02*T2+1.16000D-04*T3
      FK=1.D-03*(-35.6212+0.403494*T-1.094794D-03*T2+1.199994D-06*T3)
      RETURN
      END