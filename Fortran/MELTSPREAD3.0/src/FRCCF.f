C...................................................................... 
C.....SUBROUTINE FRCCF EVALUATES THE FRICTION COEFFICIENT GIVEN THE     
C.....REYNOLDS NUMBER AND HYDRAULIC DIAMETER.
C...................................................................... 
      SUBROUTINE FRCCF(DHYD3J,REJ,SANDDW,FRJM12)    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/PROPI/ NMATC,NMATFI,NMATFF,ICAOH2,ICACO3,IMCCO3,IFH2O,     
     1IVH2O,ICK2O,IVK2O,INA2O,ITIO2,ISIO2,ICAO,IMGO,IAL2O3,IFEO,IFE2O3, 
     2IFE3O4,IFE,ICR,INI,IZR,IU,IB4C,IB,ICR2O3,INIO,IUO2,IZRO2,IB2O3,   
     3NBCINT(999),NODES(999),NUMNOD,IMOX(28),IBAS,ILCS,ILL,NACTIV(999),     
     4NBFRZM,NBFRZO,NABLFM,NDRNFM,NTYPMT(999,999),NDOOR,NBMADJ,NUMSHV,   
     5NBSHL(999),NCRSTS(999),NUMSHH,NSKIPE,NACSH,NWAT,NPEND,NBCBOT,
     6NCRTOP,NBFZOE,NBFZME,NGEOM,NL(999,999),NFRSCT,NLOGSH,NADAB,NINVIS,
     7NSIMST,NSTEEL,NOVHT,NOVTK,NOVUM,NOVEM,NOVSIG,NPLLOC(999),NBPL,
     8IXP(999),IYP(999),NVTPE,NSOLTP,NSOLF,NINTF,IFLGA(999),ICTYPE,ICTC,
     9NPOURS,NSMP,NPED,NDOR,NSHL,NANULS,NSMPCV,NBOIL,NSWALL,NMVER,
     1NCRTEM,NTHINC,NVELP,NITMAX,NENMAX,NPRINT,NPFREQ,NTIMSPC,NPLFREQ,
     2NPLTOT,NJET,NJETP,NJETD,NJETND,NODCAP,NVELPW,NITMAXW,NENMXW,
     3IFLGJ(999),NBEDCQ,ISHELE
      IF(NINVIS.NE.1) GO TO 22
      FRJM12=0.0D0
      GO TO 40
   22 CONTINUE
      CONV=0.001D0    
      DLN10=1.D0/DLOG(10.D0)    
      REJ=DMAX1(1.D-10,REJ)
      IF(REJ.GT.2300.D0) GO TO 20  
      FRJM12=24.D0/REJ
      GO TO 40
   20 CONTINUE
      FTURBJ=0.0791D0/REJ**0.25D0
      XFRJ=1.D0/SQRT(4.D0*FTURBJ)
      RE318=18.7D0/REJ
      IFRIC=0
      DO 1219 I=1,21  
      IFRIC=IFRIC+1   
      DLNARG=2.D0*SANDDW/DHYD3J+RE318*XFRJ
      F=XFRJ-1.74D0+2.D0*DLN10*DLOG(DLNARG)
      DFDX=1.D0+(2.D0*DLN10*RE318)/DLNARG 
      DXFRJ=-F/DFDX   
      XFRJ=XFRJ+DXFRJ 
      IF(ABS(DXFRJ).LE.CONV) GO TO 30     
 1219 CONTINUE
      WRITE(2,*) 'WARNING: FRICTION FACTOR NOT FOUND, REJ=',REJ
   30 CONTINUE
      FRJM12=1.D0/(4.D0*XFRJ*XFRJ)
   40 CONTINUE
      RETURN
      END