C...................................................................... 
C.....GIVEN THE TEMPERATURE OF THE FLUID MATERIAL, THIS ROUTINE
C.....CALCULATES THE THERMAL CONDUCTIVITY OF THE OXIDE AND METAL
C.....PHASES AS WELL AS THE THERMAL CONDUCTIVITY OF THE OXIDE-
C.....METAL MIXTURE.  
C...................................................................... 
      SUBROUTINE CONDF(T,HMOL,FKO,FKM,FK) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION HMOL(28)
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
      COMMON/DATOVER/COVS,COVL,DHSOV,XMOLOV,ROVS,ROVL,TKOVS,TKOVL,VISOV,
     1EMOV,SIGOV
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
      DATA FKSOLM /2.58D0/      
      DATA FKMMLM /28.8D0/      
C.....CHECK IF INPUT DATA IS BEING USED TO OVERWRITE THE THERMAL 
C.....CONDUCTIVITY DATA AND IF SO, SET ACCORDINGLY.
      IF(NOVTK.EQ.0) GO TO 887
      IF(T.GT.TFOS) GO TO 888
      FKM=TKOVS
      FKO=TKOVS
      FK=TKOVS
      GO TO 4788 
  888 CONTINUE
      FKM=TKOVL
      FKO=TKOVL
      FK=TKOVL
      GO TO 4788
  887 CONTINUE
C.....FOR GENERAL CASE, CALCULATE THE COLLAPSED VOLUME FRACTIONS, VFLO 
C.....AND VLFM, OF THE COLLAPSED OXIDE AND METAL PHASES.
      HT=ZERO
      HTM=ZERO
      HTO=ZERO
      DO 399 I=NMATFI,NMATFF    
      HMOLI=HMOL(I)   
      HT=HT+HMOLI     
      IF(IMOX(I).NE.1) GO TO 299
      HTO=HTO+HMOLI   
      GO TO 399
  299 CONTINUE
      HTM=HTM+HMOLI   
  399 CONTINUE
      DENOM=0.D0      
      IF(HT.GT.0.D0) DENOM=ONE/HT
      VFLM=HTM*DENOM  
      VFLO=HTO*DENOM  
      IF(VFLM.GT.0.D0.OR.VFLO.GT.0.D0) GO TO 351    
      FK=FKSOLM
      GO TO 4788      
  351 CONTINUE
      TSEND=T
      IF(VFLM.LE.0.D0) GO TO 402
      CALL CONDM(TSEND,HMOL,FKM)
  402 CONTINUE
      IF(VFLO.LE.0.D0) GO TO 403
      CALL CONDO(TSEND,HMOL,FKO)
  403 CONTINUE
      IF (VFLM .GT. 0.5) GO TO 210
      VFC=VFLO
      VFD=VFLM
      FKC=FKO
      FKD=FKM
      GO TO 220
  210 VFC=VFLM
      VFD=VFLO
      FKC=FKM
      FKD=FKO
  220 CONTINUE
C.....THE FOLLOWING CALCULATES THE EFFECTIVE THERMAL CONDUCTIVITY OF THE
C.....OXIDE-METAL MIXTURE USING THE EQUATION OF WILLHITE, KUNII, AND    
C.....SMITH.
      FK=FKD
      IF(VFD.LE.0.D0) FK=FKC    
      IF(VFD.GT.0.D0.AND.VFC.GT.0.D0) FK=(FKC**VFC)*(FKD**VFD)
      FK=DMAX1(FK,FKSOLM)
      FK=DMIN1(FK,FKMMLM)
 4788 CONTINUE
      RETURN
      END