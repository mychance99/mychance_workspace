C...................................................................... 
C.....SUBROUTINE INTERP DOES A SIMPLE LINEAR INERPOLATION IN A TABLE
C.....OF X,Y DATA GIVEN THE VALUE OF X.
C...................................................................... 
      SUBROUTINE INTERP(NDIAG,NPOINT,XSEND,TC,QC,YC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TC(99),QC(99)
      DO 30 K=1,NPOINT-1
      F=(XSEND-TC(K))/(TC(K+1)-TC(K))
      IF(F.LT.0.D0) GO TO 30
      IF(F.GT.1.D0) GO TO 30
      YC=QC(K)+F*(QC(K+1)-QC(K))
      GO TO 25
   30 CONTINUE
      IF(NDIAG.EQ.1) WRITE(2,*) 'WARNING: JET DIAM. NOT FOUND IN INTERPO 
     1LATION TABLE, TIMESEND= ',XSEND
      IF(NDIAG.EQ.2) WRITE(2,*) 'WARNING: DOWNCOMER FLOW CORD LENGTH NOT
     1 FOUND IN INTERPOLATION TABLE, XSEND= ',XSEND
      IF(NDIAG.EQ.3) WRITE(2,*) 'WARNING: WATER INJECTION FLOWRATE NOT F
     1OUND IN INTERPOLATION TABLE, TIMESEND= ',XSEND
      IF(NDIAG.EQ.4) WRITE(2,*) 'WARNING: WATER INJECTION TEMPERATURE NO
     1T FOUND IN INTERPOLATION TABLE, TIMESEND= ',XSEND
   25 CONTINUE
      RETURN
      END