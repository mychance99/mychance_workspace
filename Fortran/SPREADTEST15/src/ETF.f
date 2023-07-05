C...................................................................... 
C.....GIVEN THE FLUID TEMPERATURE, THIS ROUTINE CALCULATES THE FLUID    
C.....SPECIFIC ENTHALPY AND THE DERIVATIVE OF THE TEMPERATURE WITH      
C.....RESPECT TO THE SPECIFIC ENTHALPY. THE FLUID IS ASSUMED TO
C.....CONSIST OF DISTINCT OXIDE AND METAL PHASES EACH OF WHICH IS
C.....CHARACTERIZED BY DISTINCT SOLIDUS AND LIQUIDUS TEMPERATURES.      
C...................................................................... 
      SUBROUTINE ETF(E,HIN,T,DTDE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION HIN(28),HMOL(28)
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
      COMMON/PRINTR/ QFEH2O,QCRH2O,QZRH2O,QFECO2,QCRCO2,QZRCO2,
     4XFH2OU,XFCAOH,XFMGCA,XFCACO,XZRMRE,XFEMRE,XCRMRE,XZRORE,XFEORE,   
     5XCRORE,XMH2O,XMCO2,XMCACO,XMMGCA,XMCAOH,TFWL,TFWS,TBWL,TBWS,      
     6TMCAL,TMCAS,TCAL,TCAS,TFOS,TFOL,TFMS,TFML,XVISC(28),SVISC(28),    
     7XMOL(28),FMMOL(28),ROM(28),ROMLIQ(28),AEQM(28,2),BEQM(28,2),
     8CEQM(28,2),ECL,ECS,ECAL,ECAS,EMCAL,EMCAS,EBWL,EBWS,EFWL,EFWS,    
     9STEF,GRAV,PI,TCS,TCL,CCL,CCS,RMASSL,WPCC,WPM,WPA,WPCS,ROC,RMASSS,   
     1ANGSHL,RSLAGL,RSLAGS,HNODT,VFAV,VGAV,QXAV,XWTSS(16),  
     2TSCS(2),TSCL(2),ESCS(2),ESCL(2),CCSS(2),CCSL(2),ROSTLS,ROSTLL,
     3XFRGAS,HMINC,TST(99),TSTOP(99),AINTP(99),BINTP(99),DRATIO(999),
     4XBCN(999),XDCN(999),XBLT(15),ADEC(99),BDEC(99),APOUR(16,99),
     5BPOUR(16,99),XWTC(14),BWIDTH,THCKCV,TMBOIL,TEBOIL,VFINT,ANGFAN,
     6ALPMAX,THSHL,XFCABL,XNDMIN,DVMX,DAVMX,DEAVMX,DEMX,TDC,QDCU,QDCUO2,
     7TSINJ,EINJ,DTINJ,TKINJ,ROEV,CPINJ,SURFT,VSINJ,EMINJ,TSINJO,EINJO,
     8DTINJO,TKINJO,ROEVO,CPINJO,SURFTO,VSINJO,EMINJO,XLEADE,ARSUM,
     9DBUBOX,UTRISE,TREMSH,DJET,DFALL,WEBER,FROUDE,XPSAITO,XPEPSTN,
     1FRAG,XMBEDT,XMBED(16),XLPENT,XLPENA,UJET,UFALL,HFALL,ERPV,
     2TJETT(99),DJETT(99),HWATP,XDOTJET,DVMXW,DAVMXW,DEMXW,DEAVMXW,
     3QJETW,XSTMJF,TINTSJF,ESAT,DRDOOR,DRANL,EI,QDCBUO2,QDCBU,
     4EBEDS,EBEDB
      DATA ZERO / 0.0 /
      DATA ONE / 1.0 /
      DATA TWO / 2.0 /
      NMATFI=8
      NMATFF=28
C.....DEFINE THE COLLAPSED DEPTHS OF THE FLUID MATERIALS HERE GIVEN THE 
C.....COLLAPSED DEPTHS, HIN, CALCULATED BY THE FLUID DYNAMICS 
C.....ROUTINES.
      DNOM=ZERO
      DO 399 I=NMATFI,NMATFF    
      HMOL(I)=HIN(I)*ROMLIQ(I)/FMMOL(I)   
      DNOM=DNOM+HMOL(I)
  399 CONTINUE
      RDNOM=ZERO      
      IF (DNOM .NE. ZERO) RDNOM=ONE/DNOM  
      DO 499 I=NMATFI,NMATFF    
      HMOL(I)=HMOL(I)*RDNOM     
  499 CONTINUE
C.....TEMPERATURE TABLE
C.....TFOL:  FLUID OXIDE PHASE LIQUIDUS   
C.....TFOS:  FLUID OXIDE PHASE SOLIDUS    
C.....TFML:  FLUID METAL PHASE LIQUIDUS   
C.....TFMS:  FLUID METAL PHASE SOLIDUS    
C.....CALCULATE THE SPECIFIC ENTHALPY AND DERIVATIVE OF THE TEMPERATURE 
C.....WITH RESPECT TO THE SPECIFIC ENTHALPY CORRESPONDING TO THE
C.....CURRENT FLUID TEMPERATURE.
      A=ZERO
      B=ZERO
      C=ZERO
      DNOM=ZERO
      DO 999 I=NMATFI,NMATFF    
      HMOLI=HMOL(I)   
      DNOM=DNOM+HMOLI*FMMOL(I)  
      J=1
      IF(IMOX(I).EQ.0) GO TO 905
      IF(T.LE.TFOS) GO TO 910
      J=2
      IF(T.GE.TFOL) GO TO 910
      ELI=CEQM(I,2)+TFOL*(BEQM(I,2)+TFOL*AEQM(I,2)) 
      ESI=CEQM(I,1)+TFOS*(BEQM(I,1)+TFOS*AEQM(I,1)) 
      B=B+((ELI-ESI)/(TFOL-TFOS))*HMOLI   
      C=C+(ESI-TFOS*((ELI-ESI)/(TFOL-TFOS)))*HMOLI
      GO TO 999
  905 CONTINUE
      IF(T.LE.TFMS) GO TO 910
      J=2
      IF(T.GE.TFML) GO TO 910
      ELI=CEQM(I,2)+TFML*(BEQM(I,2)+TFML*AEQM(I,2)) 
      ESI=CEQM(I,1)+TFMS*(BEQM(I,1)+TFMS*AEQM(I,1)) 
      B=B+((ELI-ESI)/(TFML-TFMS))*HMOLI   
      C=C+(ESI-TFMS*((ELI-ESI)/(TFML-TFMS)))*HMOLI
      GO TO 999
  910 CONTINUE     
      A=A+HMOLI*AEQM(I,J)
      B=B+HMOLI*BEQM(I,J)
      C=C+HMOLI*CEQM(I,J)
  999 CONTINUE
      RDNOM=ZERO      
      IF (DNOM .NE. ZERO) RDNOM=ONE/DNOM  
      A=A*RDNOM
      B=B*RDNOM
      C=C*RDNOM
      IF (ABS(A) .GT. 1.D-05) GO TO 1100
      E=B*T+C
      DEDT=ABS(B)     
      DTDE=ABS(ONE/DEDT)
      GO TO 1200      
 1100 E=C+T*(B+T*A)   
      DEDT=ABS(B+A*TWO*T)
      DTDE=ABS(ONE/DEDT)
 1200 CONTINUE
      RETURN
      END