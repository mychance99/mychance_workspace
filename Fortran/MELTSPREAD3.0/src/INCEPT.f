C...................................................................... 
C.....SUBROUTINE INCEPT CHECKS MELT CONSTITUENT INVENTORY AND MELT-
C.....BASEMAT INTERFACIAL CONDITIONS FOR INCEPTION OF CRUSTING, 
C.....ABLATION, OR SIMULTANEOUS CRUSTING AND ABLATION.    
C...................................................................... 
      SUBROUTINE INCEPT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION HSEND(28)
      COMMON/CATCHR/A(999),B(999),C(999),D(999),P(999),Q(999),
     1RECON(999),DELVEL(999),DSRDT(999),DELHTS(999),EMAX(999,999),
     2DBLKEN(999),DELHIT(16,999),SOURCE(16,999),SRCTOT(999),
     3HCROLD(16,999),DCDOT(999),HCRUST(16,999),DFILM(16,999),
     4DFMOLD(16,999),DENOLD(999),SIGOX(999),SIGCON(999),BET(999),
     5TFZX(999),TFZC(999),TKX(999),TKM(999),THETO(999,999),
     6THET1(999,999),OMEG0(4,999),OMEG1(4,999),TAO0(999),TAO1(999),
     7SIG0(999),SIG1(999),TARGB(999),DADOT(999),DADOTG(999),
     8XLOLD(4,999),VISREX(999),TATMS(999),HTMP(999),DHTMP(999),
     9HTMOLD(999)
      COMMON/CATCHI/NUMS(4,999),NUMOLD(4,999),NTRACK(16),NSUMP(999),
     1NPASSG(999),NCRSTM(999),NCRRT(999),NPASS(999),NBURN,NBURNO     
      COMMON/PRINTB/ELEVAT(999),ELOLD(999),XDIST(999,999),XBTW(999,999),
     1TEMP(999,999),ENTHP(999,999),EOLD(999,999),ENBLK(999),EBKOLD(999),
     2TBULK(999),HITE(16,999),HITOLD(16,999),HTOT(999),HTOLD(999),
     3HTCFT(999),HTCOEF(999),QFLUXT(999),QFLUXB(999),VEL(999),ELO(999),
     4VELOLD(999),SRSCOR(16),AREA(999),RAD(999),ARC(999),VOLCN(16),
     5VOID(999),VGJ(999),ZABLAT(999),ZABOLD(999),DCRUST(999),
     6DCROLD(999),DABCON(999),DABOLD(999),DFILMT(999),DFOLT(999),
     7SMFLX(4,999),XLSMF(4,999),TOTVOL,XFACMS(999),XMFLXA,XMCORT,
     8VCORT,TOTOX,TOTMET,QFLXT,QFLXB,TIME,DTIME,XMCOR(16),VCOR(16),
     9TCONI,RSAND,HDOWNC,TBOUND,EMISCN,PDRYWL,XDISTO(999),QOXT(999),
     1XLENSH,XBTWO(999),RCOMP,WDOOR,RSUMP,RSHELL,RPED,TPED,ELSMP,
     1TENDP(10,999),TFRZSH,TDEBRS,TKDEBR,PDEBR,CPDEBR,ENENDP(10,999),
     2ENOLDP(10,999),DXVERT,DXSNK,XDSTE(10,999),XBTE(10,999),HXLA(999),
     3HXLB(999),TSFEB(999),DCRS(999),DCRSLD(999),HCRS(16,999),
     4HCRSLD(16,999),THETE0(10,999),THETE1(10,999),QSHELL,QSHELE,TSHELI,
     5SIGOXE(999),FKOXE(999),RINJC,XLSEC,ANGSEC,ANINJC,XLCHAN,WCHNL,
     6TEFZX(999),TIMSPC(999),DXNODE(999),TNORM(999),FRCSOL(999),
     7ALPSPR(999),CRAMCON,HINTF,TSHLMX,XFRMX(999),XFROX(999),XFRTX(999),
     8XMLMX(999),XMLOX(999),XMLTX(999),XTOTX(999),TIMEO,TMAX,EDOWN(999),
     9HCP(999)
      COMMON/HCONS/HMETAL(999),HOXIDE(999)  
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
      COMMON/BLKHT/TKMLT,PM,CM,UM,EMISM,SIGMAM,TKNC,PNC,CNC,UNC    
      COMMON/BASPR/ HBMT(999,999,16),HBMOLD(999,999,16) 
      COMMON/PROPB/CPBM(999),FKB(999,999),DENB(999,999),DTDEB(999,999),
     1EMIBM(999),DENCRS(999),ENCRS(999),FKRF(999),DENRF(999),CPRF(999) 
C.....RUN THROUGH THE MESH AND UPDATE BOUNDARY CONDITIONS....
      NCHNG=0
      DO 7750 KN=1,NUMNOD
      IF(NACTIV(KN).EQ.0.OR.NADAB.EQ.1) GO TO 7750      
      IF(HTOT(KN).LE.0.D0) GO TO 7750     
      IF(NGEOM.EQ.2) GO TO 1919
      IF(KN.GT.NSMP) GO TO 1919
      IF(NSMPCV.EQ.0) GO TO 1919
      IF(NBOIL.EQ.0) GO TO 1919
      IF(NBURNO.EQ.0) GO TO 7750
 1919 CONTINUE
      NCHNG=1
      IF(NBFRZO.EQ.1.AND.NBFRZM.EQ.1) GO TO 7751    
C.....CHECK MELT CONSTITUENT INVENTORY FOR THE POSSIBILITY OF STARTING  
C.....TO CRUST EITHER METAL OR OXIDE CRUST.
      IF(NBCINT(KN).EQ.1) GO TO 7769      
      IF(NBCINT(KN).EQ.4) GO TO 7769      
      IF(NBCINT(KN).EQ.5) GO TO 7769      
      GO TO 7751      
 7769 CONTINUE
      HFROX=HOXIDE(KN)/HTOT(KN) 
      HFRM=HMETAL(KN)/HTOT(KN)  
      NCRWTC=0
      IF(NBFRZO.EQ.1) GO TO 1884
      IF(HFROX.LT.0.5D0) GO TO 1884
      IF(TBULK(KN).GT.TFOS.AND.TEMP(1,KN).LT.TFOS) NCRWTC=1   
      IF(NCRWTC.EQ.1) GO TO 7772
 1884 CONTINUE
      IF(NBFRZM.EQ.1) GO TO 7751
      IF(HMETAL(KN).GT.0.D0.AND.TEMP(1,KN).LT.TFMS) NCRWTC=2  
      IF(NCRWTC.EQ.0) GO TO 7751
 7772 CONTINUE
C.....MELT CONSTITUENCY IS SUFFICIENT FOR INCEPTION OF CRUST GROWTH.    
C.....CHECK THERMAL CONDITIONS. 
      GO TO(7773,7774),NCRWTC   
 7773 CONTINUE
      TFREEZ=TFOS     
      TFRZL=TFOL      
      HCOMP=HOXIDE(KN)
      GO TO 7775      
 7774 CONTINUE
      TFREEZ=TFMS     
      TFRZL=TFML      
      HCOMP=HMETAL(KN)
 7775 CONTINUE
      IF(NBCINT(KN).EQ.5) GO TO 7776      
      IF(TEMP(1,KN).LT.TFREEZ) GO TO 7777 
      IF(NBCINT(KN).EQ.4) GO TO 7750      
      GO TO 7751      
 7776 CONTINUE
      ZCOMP=XDISTO(KN)-ZABLAT(KN)
      IP=2  
      IF(ZCOMP.LT.XBTWO(KN)) IP=1
      IF(NTYPMT(IP,KN).GT.3) GO TO 7780   
      TMLTL=TSCL(NL(IP,KN))     
      GO TO 7784      
 7780 CONTINUE
      IROUT=NTYPMT(IP,KN)-3     
      GO TO(7782,7783),IROUT    
 7782 CONTINUE
      TMLTL=TFML      
      GO TO 7784      
 7783 CONTINUE
      TMLTL=TFOL      
 7784 CONTINUE
      CALL ASINEA(DFILM,KN,HSEND,ROMLIQ,WPCC,WPM,WPA,WPCS,ROC,RSLAGL,0) 
      CALL CONDF(TMLTL,HSEND,FKO,FKM,TKMLT)
      C1=HTCOEF(KN)+TKMLT/DFILMT(KN)      
      C2=HTCOEF(KN)*TBULK(KN)+(TKMLT*TEMP(1,KN))/DFILMT(KN)   
      TINT=C2/C1      
      IF(TINT.GE.TFREEZ) GO TO 7750
 7777 CONTINUE
C.....CRUST GROWTH CAN BEGIN.  FIND THE THERMAL PROPERTIES OF THE
C.....CRUSTING MATERIAL AND BASEMAT.      
      NCRSTM(KN)=NCRWTC
      IF(NTHINC.NE.1) GO TO 8111
      IF(NBCINT(KN).GT.1) GO TO 7788      
      NBCINT(KN)=3    
      GO TO 7750      
 8111 CONTINUE
      CALL ASINEA(HITE,KN,HSEND,ROMLIQ,WPCC,WPM,WPA,WPCS,ROC,RSLAGL,NCRS
     1TM(KN))
      CALL ETF(EFL,HSEND,TFRZL,DTDE)      
      CALL ASINEA(HITE,KN,HSEND,ROM,WPCC,WPM,WPA,WPCS,ROC,RSLAGS,NCRSTM(
     1KN))  
      CALL ETF(EFS,HSEND,TFREEZ,DTDE)     
      HOX=EFL-EFS     
      TEV=TFREEZ-1.D0 
      CALL ETF(EPS,HSEND,TEV,DTDE)
      COX=1.D0/DTDE   
      CALL CONDF(TFREEZ,HSEND,FKO,FKM,TKOX)
      CALL MASFRA(HITE,ROM,KN,RSLAGS,ROC,DENOX,NCRSTM(KN))    
      CALL MASFRA(HITE,ROMLIQ,KN,RSLAGL,ROC,DENOXL,NCRSTM(KN))
      IF(NBCINT(KN).GT.1) GO TO 7788      
      IF(NTYPMT(1,KN).GT.3) GO TO 7789    
      TFRZS=TSCS(NL(1,KN))      
      ESOL=ESCS(NL(1,KN))
      DTDESL=1.D0/CCSS(NL(1,KN))
      HABL=ESCL(NL(1,KN))-ESCS(NL(1,KN))  
      DRATLS=RSLAGL/RMASSS      
      IF(NL(1,KN).EQ.2) DRATLS=ROSTLL/ROSTLS
      DENML=ROC*RSLAGL
      IF(NL(1,KN).EQ.2) DENML=ROSTLL     
      CPML=CCSL(NL(1,KN))
      CALL CONDC(TSCL(NL(1,KN)),NTYPMT(1,KN),TKML,NL(1,KN))   
      CALL CONDC(TEMP(1,KN),NTYPMT(1,KN),TK1,NL(1,KN))
      CALL TEC(ENTHP(1,KN),EPS,TPS,DTDE1,RMASS1,YCACO3,YMCCO3,YCAOH2,YFH
     12O,NL(1,KN))    
      CP1=1.D0/DTDE1  
      DEN1=ROC
      IF(NL(1,KN).EQ.2) DEN1=ROSTLS      
      GO TO 7790      
 7789 CONTINUE
      IROUT=NTYPMT(1,KN)-3      
      GO TO(6333,6533),IROUT    
 6333 CONTINUE
      TFRZL=TFML      
      TFRZS=TFMS      
      GO TO 6534      
 6533 CONTINUE
      TFRZL=TFOL      
      TFRZS=TFOS      
 6534 CONTINUE
      CALL ASINED(1,KN,HSEND,ROMLIQ,WPCC,WPM,WPA,WPCS,ROC,RSLAGL,0)     
      CALL CONDF(TFRZL,HSEND,FKO,FKM,TKML)
      CALL ETF(ELIQ,HSEND,TFRZL,DTDE1)    
      TEV=TFRZL-1.D0  
      CALL ETF(EPS,HSEND,TEV,DTDE)
      CPML=1.D0/DTDE  
      CALL MASFRD(ROMLIQ,1,KN,RSLAGL,ROC,DENML)     
      CALL ASINED(1,KN,HSEND,ROM,WPCC,WPM,WPA,WPCS,ROC,RSLAGS,0)
      CALL ETF(ESOL,HSEND,TFRZS,DTDE)     
      HABL=ELIQ-ESOL  
      DTDESL=DTDE     
      CALL MASFRD(ROM,1,KN,RSLAGS,ROC,DENS)
      DRATLS=DENML/DENS
      CALL ASINED(1,KN,HSEND,ROM,WPCC,WPM,WPA,WPCS,ROC,RSLAGS,0)
      CALL CONDF(TEMP(1,KN),HSEND,FKO,FKM,TK1)      
      CALL ETF(EPASS,HSEND,TEMP(1,KN),DTDE1)
      CP1=1.D0/DTDE1  
      CALL MASFRD(ROM,1,KN,RSLAGS,ROC,DEN1)
 7790 CONTINUE
      IF(NTYPMT(2,KN).GT.3) GO TO 7791    
      CALL CONDC(TEMP(2,KN),NTYPMT(2,KN),TK2,NL(2,KN))
      GO TO 7792      
 7791 CONTINUE
      CALL ASINED(2,KN,HSEND,ROM,WPCC,WPM,WPA,WPCS,ROC,RSLAGS,0)
      CALL CONDF(TEMP(2,KN),HSEND,FKO,FKM,TK2)      
 7792 CONTINUE
      TKE=(XDIST(1,KN)*TK1*TK2)/(TK2*XBTW(1,KN)+TK1*(XDIST(1,KN)-XBTW(1,
     1KN))) 
C.....CALCULATE INITIAL CRUST DEPTH WHEN CRUST GROWTH STARTS FROM THE   
C.....CONDITION NBC=1 (NO CRUSTING, NO ABLATION).   
      ALPHA=(TKE*DENOX*HOX)/(XDIST(1,KN)*CP1*DEN1*XBTW(1,KN)) 
      BETA=TKOX*(TFREEZ-TEMP(1,KN))
      GAMMA=HTCOEF(KN)*(TBULK(KN)-TEMP(1,KN))+(TKOX*DENOX*HOX)/(CP1*DEN1
     1*XBTW(1,KN))    
      Z1=(GAMMA/ALPHA)*(GAMMA/ALPHA)+(4.D0*BETA)/ALPHA
      DCRUST(KN)=-0.5D0*(GAMMA/ALPHA)+0.5D0*SQRT(Z1)
C.....CHECK TO DETERMINE IF THE INITIAL CRUST MASS EXCEEDS THAT
C.....AVAILABLE IN THE MELT LAYER.  IF SO, MAKE APPROPRIATE ADJUSTMENT. 
      HAVAIL=(DENOXL*HCOMP)/DENOX
      IF(DCRUST(KN).LE.HAVAIL) GO TO 7793 
      DCRUST(KN)=HAVAIL
 7793 CONTINUE
C.....DEFICIT THE BULK MELT THE AMOUNT OF THE INITIAL CRUST MASS, AND   
C.....ADJUST BASEMAT ELEVATION. 
      IF(NCRSTM(KN).EQ.2) GO TO 7794     
      KINL=1
      KINH=7
      KLOW=8
      KHIGH=16
      GO TO 7795      
 7794 CONTINUE
      KINL=8
      KINH=16
      KLOW=1
      KHIGH=7
 7795 CONTINUE
      DO 7796 KTP=KLOW,KHIGH    
      CALL INDEX(KTP,KJ)
      HCRUST(KTP,KN)=(HITE(KTP,KN)*DENOX*DCRUST(KN)*ROMLIQ(KJ))/(HCOMP*D
     1ENOXL*ROM(KJ))  
      HITE(KTP,KN)=HITE(KTP,KN)*(1.D0-(DCRUST(KN)*DENOX)/(DENOXL*HCOMP))
      HITE(KTP,KN)=DMAX1(HITE(KTP,KN),0.D0)
 7796 CONTINUE
      DO 6371 KTP=KINL,KINH     
      HCRUST(KTP,KN)=0.D0
 6371 CONTINUE
      HTOT(KN)=0.D0   
      HMETAL(KN)=0.D0 
      HOXIDE(KN)=0.D0 
      DO 7797 KTP=1,16
      HTOT(KN)=HTOT(KN)+HITE(KTP,KN)      
      IF(KTP.LE.7) HMETAL(KN)=HMETAL(KN)+HITE(KTP,KN)
      IF(KTP.GE.8) HOXIDE(KN)=HOXIDE(KN)+HITE(KTP,KN)
 7797 CONTINUE
      ELEVAT(KN)=ELEVAT(KN)+DCRUST(KN)    
C.....CALCULATE NEW CRUST/BASEMAT INTERFACIAL TEMPERATURE.  IF IT
C.....EXCEEDS BASEMAT MELTING TEMPERATURE, START ABLATION CALCULATION.  
      IF(NPASS(KN).EQ.1) GO TO 7750
      ECOMP=ENTHP(1,KN)+(DENOX*HOX*DCRUST(KN))/(DEN1*XBTW(1,KN))
      IF(ECOMP.GT.ESOL) GO TO 3621
      NBCINT(KN)=2    
      ENTHP(1,KN)=ECOMP
      IF(NTYPMT(1,KN).GT.3) GO TO 7798    
      CALL TEC(ENTHP(1,KN),EMAX(1,KN),TEMP(1,KN),DTDE,RMASS,YCACO3,YMCCO
     13,YCAOH2,YFH2O,NL(1,KN))  
      DTDEB(1,KN)=DTDE
      GO TO 7750      
 7798 CONTINUE
      CALL ASINED(1,KN,HSEND,ROM,WPCC,WPM,WPA,WPCS,ROC,RSLAGS,0)
      CALL ETF(ENTHP(1,KN),HSEND,TEMP(1,KN),DTDE)   
      DTDEB(1,KN)=DTDE
      GO TO 7750      
 3621 CONTINUE
C.....CRUSTING WITH ABLATION STARTS.  INITIALIZE ABLATION CALCULATION.  
      NBCINT(KN)=6    
      IF(NDRNFM.EQ.2) NBCINT(KN)=8
      IF(NTHINC.EQ.1.AND.NBCINT(KN).EQ.6) NBCINT(KN)=7
      IF(NTHINC.EQ.1.AND.NBCINT(KN).EQ.8) NBCINT(KN)=9
      DFS=(XBTW(1,KN)*(ECOMP-ESOL))/HABL  
      DFSMAX=XBTW(1,KN)*(1.D0-XFCABL)     
      DFS=DMIN1(DFS,DFSMAX)     
      DFM=DFS/DRATLS  
      IF(NBCINT(KN).EQ.8.OR.NBCINT(KN).EQ.9) GO TO 7799
C.....POROUS CRUST ABLATION CASE TREATED HERE.  ADD INITIALLY ABLATED   
C.....BASEMAT MATERIAL TO BULK MELT.      
      HTOT(KN)=HTOT(KN)+DFM     
      IF(NTYPMT(1,KN).GT.3) GO TO 7900    
      IF(NL(1,KN).EQ.2) GO TO 7900
      HITE(16,KN)=HITE(16,KN)+DFM
      HOXIDE(KN)=HOXIDE(KN)+DFM 
      GO TO 7901      
 7900 CONTINUE
      DO 7902 KTP=1,16
      CALL INDEX(KTP,KJ)
      DENR=ROM(KJ)/ROMLIQ(KJ)   
      IF(KTP.EQ.16) DENR=RMASSS/RSLAGL    
      DINCR=(DFS*HBMT(1,KN,KTP))/XBTW(1,KN)
      IF(NL(1,KN).EQ.2) DINCR=DFS*XWTSS(KTP)
      HITE(KTP,KN)=HITE(KTP,KN)+DINCR*DENR
      HBMT(1,KN,KTP)=HBMT(1,KN,KTP)-DINCR 
      HBMT(1,KN,KTP)=DMAX1(HBMT(1,KN,KTP),0.D0)     
 7902 CONTINUE
      HMETAL(KN)=0.D0 
      HOXIDE(KN)=0.D0 
      DO 6077 KTP=1,16
      IF(KTP.LE.7) HMETAL(KN)=HMETAL(KN)+HITE(KTP,KN)
      IF(KTP.GE.8) HOXIDE(KN)=HOXIDE(KN)+HITE(KTP,KN)
 6077 CONTINUE
      GO TO 7901      
 7799 CONTINUE
C.....NON-POROUS CRUST CASE TREATED HERE.  INITIALIZE MOLTEN FILM
C.....INVENTORY.      
      DFILMT(KN)=DFM  
      IF(NTYPMT(1,KN).GT.3) GO TO 8010    
      IF(NL(1,KN).EQ.2) GO TO 8010
      DO 8011 KTP=1,16
      DFILM(KTP,KN)=0.D0
      IF(KTP.LT.16) GO TO 8011  
      DFILM(KTP,KN)=DFM
 8011 CONTINUE
      GO TO 7901      
 8010 CONTINUE
      DO 8012 KTP=1,16
      CALL INDEX(KTP,KJ)
      DENR=ROM(KJ)/ROMLIQ(KJ)   
      IF(KTP.EQ.16) DENR=RSLAGS/RSLAGL    
      DINCR=(DFS*HBMT(1,KN,KTP))/XBTW(1,KN)
      IF(NL(1,KN).EQ.2) DINCR=DFS*XWTSS(KTP)
      DFILM(KTP,KN)=DENR*DINCR  
      IF(NL(1,KN).EQ.2) GO TO 8012
      HBMT(1,KN,KTP)=HBMT(1,KN,KTP)-DINCR 
      HBMT(1,KN,KTP)=DMAX1(HBMT(1,KN,KTP),0.D0)     
 8012 CONTINUE
 7901 CONTINUE
      DTDEB(1,KN)=DTDESL
      TEMP(1,KN)=TFRZS
      ENTHP(1,KN)=ESOL
      EMAX(1,KN)=ENTHP(1,KN)    
      ZABLAT(KN)=XDIST(1,KN)-DFS
      ZABOLD(KN)=ZABLAT(KN)
      XDISTO(KN)=XDIST(1,KN)    
      XBTWO(KN)=XBTW(1,KN)      
      ELEVAT(KN)=ELEVAT(KN)-DFS 
      IF(NBCINT(KN).EQ.8.OR.NBCINT(KN).EQ.9) ELEVAT(KN)=ELEVAT(KN)+DFM  
      IF(NTYPMT(1,KN).LE.3) DABCON(KN)=DABOLD(KN)+DFS
      ELOLD(KN)=ELEVAT(KN)
      DABOLD(KN)=DABCON(KN)
      GO TO 7750      
 7788 CONTINUE
C.....INITIALIZE CRUST GROWTH CALCULATION WHEN CRUST GROWTH STARTS OVER 
C.....ABLATING BASEMAT.  FIRST FIND FILM & BASEMAT THERMAL PROPERTIES.  
      ZCOMP=XDISTO(KN)-ZABLAT(KN)
      IF(ZCOMP.LT.XBTWO(KN)) GO TO 3778   
      ZLENGT=ZABLAT(KN)+XBTW(2,KN)
      IP=2  
      GO TO 3779      
 3778 CONTINUE
      ZLENGT=XBTWO(KN)-(XDISTO(KN)-ZABLAT(KN))      
      IP=1  
 3779 CONTINUE
C.....SPECIFY ABLATION (MELTING) TEMPERATURE & MOLTEN FILM PROPERTIES.  
      IF(NTYPMT(IP,KN).GT.3) GO TO 3780   
      TFRZBM=TSCS(NL(IP,KN))    
      ESOL=ESCS(NL(IP,KN))      
      HABL=ESCL(NL(IP,KN))-ESCS(NL(IP,KN))
      DENS=ROC
      IF(NL(IP,KN).EQ.2) DENS=ROSTLS     
      DRATLS=RSLAGL/RMASSS      
      IF(NL(IP,KN).EQ.2) DRATLS=ROSTLL/ROSTLS      
      IF(NBCINT(KN).EQ.5) GO TO 3990      
      DENF=ROC*RSLAGL 
      IF(NL(IP,KN).EQ.2) DENF=ROSTLL     
      CPF=CCSL(NL(IP,KN))
      CALL CONDC(TSCL(NL(IP,KN)),NTYPMT(IP,KN),TKF,NL(IP,KN)) 
      GO TO 3781      
 3780 CONTINUE
      IF(NTYPMT(IP,KN).EQ.5) GO TO 3782   
      TFRZBM=TFMS     
      TFRZBL=TFML     
      GO TO 3783      
 3782 CONTINUE
      TFRZBM=TFOS     
      TFRZBL=TFOL     
 3783 CONTINUE
      CALL MASFRD(ROMLIQ,IP,KN,RSLAGL,ROC,DENL)     
      CALL MASFRD(ROM,IP,KN,RSLAGS,ROC,DENS)
      DRATLS=DENL/DENS
      CALL ASINED(IP,KN,HSEND,ROM,WPCC,WPM,WPA,WPCS,ROC,RSLAGS,0)
      CALL ETF(ESOL,HSEND,TFRZBM,DTDE)    
      CALL ASINED(IP,KN,HSEND,ROMLIQ,WPCC,WPM,WPA,WPCS,ROC,RSLAGL,0)    
      CALL ETF(ELIQ,HSEND,TFRZBL,DTDEF)   
      HABL=ELIQ-ESOL  
      IF(NBCINT(KN).EQ.5) GO TO 3990      
      DENF=DENL
      TEV=TFRZBL+1.D0 
      CALL ETF(EPS,HSEND,TEV,DTDEO)
      CPF=1.D0/DTDEO  
      CALL CONDF(TFRZBL,HSEND,FKO,FKM,TKF)
      GO TO 3781      
 3990 CONTINUE
      CALL MASFRA(DFILM,ROMLIQ,KN,RSLAGL,ROC,DENF,0)
      CALL ASINEA(DFILM,KN,HSEND,ROMLIQ,WPCC,WPM,WPA,WPCS,ROC,RSLAGL,0) 
      CALL ETF(EPASS,HSEND,TFRZBL,DTDE)   
      CPF=1.D0/DTDE   
      CALL CONDF(TFRZBL,HSEND,FKO,FKM,TKF)
 3781 CONTINUE
C.....CALCULATE INITIAL CRUST DEPTH FOR CASE IN WHICH STATIONARY MOLTEN 
C.....FILM PRE-EXISTS ON THE BASEMAT.     
      IF(NBCINT(KN).EQ.4) GO TO 3791      
      IF(NDRNFM.EQ.1) GO TO 3792
      NBCINT(KN)=8    
      IF(NTHINC.EQ.1) NBCINT(KN)=9
      GO TO 3793      
 3792 CONTINUE
      NBCINT(KN)=6    
      IF(NTHINC.EQ.1) NBCINT(KN)=7
 3793 CONTINUE
      IF(NTHINC.EQ.1) GO TO 7677
      ALPHA=(2.D0*TKF*DENOX*HOX)/(CPF*DENF*DFILMT(KN)*DFILMT(KN))
      CON1=HTCOEF(KN)*TBULK(KN)+(TKF*TFRZBM)/DFILMT(KN)
      CON2=HTCOEF(KN)+TKF/DFILMT(KN)      
      TINT=CON1/CON2  
      BETA=(2.D0*DENOX*HOX*TKOX)/(CPF*DENF*DFILMT(KN))+(TKF*(TINT-TFRZBM
     1))/DFILMT(KN)   
      GAMMA=TKOX*(TFREEZ-TINT)  
      CON3=(BETA/ALPHA)*(BETA/ALPHA)+(4.D0*GAMMA)/ALPHA
      DCRUST(KN)=-(0.5D0*BETA)/ALPHA+0.5D0*SQRT(CON3)
C.....CHECK TO SEE IF THIS AMOUNT OF CRUST GROWTH EXCEEDS AVAILABLE IN  
C.....MELT LAYER.  IF SO, MAKE APPROPRIATE ADJUSTMENT.
      HCOMPS=(HCOMP*DENOXL)/DENOX
      IF(DCRUST(KN).LE.HCOMPS) GO TO 3794 
      DCRUST(KN)=HCOMPS
 3794 CONTINUE
C.....INITIALIZE CRUST MASS INVENTORY & DEFICIT BULK MELT THE AMOUNT OF 
C.....CRUSTED MATERIAL.
      IF(NCRSTM(KN).EQ.2) GO TO 3795     
      KINL=1
      KINH=7
      KLOW=8
      KHIGH=16
      GO TO 3796      
 3795 CONTINUE
      KINL=8
      KINH=16
      KLOW=1
      KHIGH=7
 3796 CONTINUE
      DO 3797 KTP=KLOW,KHIGH    
      CALL INDEX(KTP,KJ)
      HCRUST(KTP,KN)=(HITE(KTP,KN)*DENOX*DCRUST(KN)*ROMLIQ(KJ))/(HCOMP*D
     1ENOXL*ROM(KJ))  
      HITE(KTP,KN)=HITE(KTP,KN)*(1.D0-(DCRUST(KN)*DENOX)/(DENOXL*HCOMP))
      HITE(KTP,KN)=DMAX1(HITE(KTP,KN),0.D0)
 3797 CONTINUE
      DO 3798 KTP=KINL,KINH     
      HCRUST(KTP,KN)=0.D0
 3798 CONTINUE
      HMETAL(KN)=0.D0 
      HOXIDE(KN)=0.D0 
      HTOT(KN)=0.D0   
      DO 7691 KTP=1,16
      IF(KTP.LE.7) HMETAL(KN)=HMETAL(KN)+HITE(KTP,KN)
      IF(KTP.GE.8) HOXIDE(KN)=HOXIDE(KN)+HITE(KTP,KN)
      HTOT(KN)=HTOT(KN)+HITE(KTP,KN)      
 7691 CONTINUE
 7677 CONTINUE
C.....ADJUST ELEVATION TO ACCOUNT FOR CRUSTED MATERIAL.
      IF(NBCINT(KN).EQ.8.OR.NBCINT(KN).EQ.9) GO TO 3830
      ELEVAT(KN)=ELEVAT(KN)+DCRUST(KN)-DFILMT(KN)   
      GO TO 3831      
 3830 CONTINUE
      ELEVAT(KN)=ELEVAT(KN)+DCRUST(KN)    
 3831 CONTINUE
C.....DRAIN THE FILM INTO BULK MELT IF NEW CRUST IS POROUS.   
      IF(NBCINT(KN).EQ.8.OR.NBCINT(KN).EQ.9) GO TO 7750
      DFILMT(KN)=0.D0 
      HTOT(KN)=0.D0   
      HMETAL(KN)=0.D0 
      HOXIDE(KN)=0.D0 
      DO 3832 KTP=1,16
      IF(KTP.LE.7) HMETAL(KN)=HMETAL(KN)+DFILM(KTP,KN)
      IF(KTP.GE.8) HOXIDE(KN)=HOXIDE(KN)+DFILM(KTP,KN)
      HITE(KTP,KN)=HITE(KTP,KN)+DFILM(KTP,KN)
      HTOT(KN)=HTOT(KN)+HITE(KTP,KN)      
      DFILM(KTP,KN)=0.D0
 3832 CONTINUE
      GO TO 7750      
 3791 CONTINUE
C.....CALCULATE INITIAL CRUST DEPTH FOR CASES IN WHICH MOLTEN FILM      
C.....WAS PREVIOUSLY CONTINUOUSLY MIXED OFF SURFACE.  FIRST SET NEW     
C.....BOUNDARY CONDITION ACCORDING TO WHETHER OR NOT CRUST IS POROUS.   
      IF(NDRNFM.EQ.1) GO TO 3833
      NBCINT(KN)=8    
      IF(NTHINC.EQ.1) NBCINT(KN)=9
      GO TO 3834      
 3833 CONTINUE
      NBCINT(KN)=6    
      IF(NTHINC.EQ.1) NBCINT(KN)=7
 3834 CONTINUE
      IF(NTHINC.EQ.1) GO TO 3417
C.....CALCULATE CRUST DEPTH BASED ON CONVECTIVE EXCHANGE ARGUMENT.      
      DFS=DTIME*ABS(DADOTG(KN)) 
      DCRUST(KN)=(DFS*DENS*CPF*(TFREEZ-TFRZBM))/(DENOX*HOX)   
C.....CHECK TO SEE IF CRUST GROWTH EXCEEDS THAT AVAILABLE IN MELT LAYER.
      HCOMPS=(HCOMP*DENOXL)/DENOX
      IF(DCRUST(KN).LE.HCOMPS) GO TO 3835 
      DCRUST(KN)=HCOMPS
 3835 CONTINUE
C.....INITIALIZE CRUST MASS INVENTORY & DEFICIT BULK MELT INVENTORY.    
      IF(NCRSTM(KN).EQ.2) GO TO 3836     
      KINL=1
      KINH=7
      KLOW=8
      KHIGH=16
      GO TO 3837      
 3836 CONTINUE
      KINL=8
      KINH=16
      KLOW=1
      KHIGH=7
 3837 CONTINUE
      DO 3838 KTP=KLOW,KHIGH    
      CALL INDEX(KTP,KJ)
      HCRUST(KTP,KN)=(HITE(KTP,KN)*DENOX*DCRUST(KN)*ROMLIQ(KJ))/(HCOMP*D
     1ENOXL*ROM(KJ))  
      HITE(KTP,KN)=HITE(KTP,KN)*(1.D0-(DCRUST(KN)*DENOX)/(DENOXL*HCOMP))
      HITE(KTP,KN)=DMAX1(HITE(KTP,KN),0.D0)
 3838 CONTINUE
      DO 3839 KTP=KINL,KINH     
      HCRUST(KTP,KN)=0.D0
 3839 CONTINUE
      HTOT(KN)=0.D0   
      HMETAL(KN)=0.D0 
      HOXIDE(KN)=0.D0 
      DO 3799 KTP=1,16
      HTOT(KN)=HTOT(KN)+HITE(KTP,KN)      
      IF(KTP.LE.7) HMETAL(KN)=HMETAL(KN)+HITE(KTP,KN)
      IF(KTP.GE.8) HOXIDE(KN)=HOXIDE(KN)+HITE(KTP,KN)
 3799 CONTINUE
      IF(NBCINT(KN).EQ.8.OR.NBCINT(KN).EQ.9) GO TO 3417
      ELEVAT(KN)=ELEVAT(KN)+DCRUST(KN)    
      GO TO 7750      
 3417 CONTINUE
C.....INITIALIZE MOLTEN FILM DEPTH & INVENTORY FOR CASE OF NON-POROUS   
C.....CRUST GROWTH.   
      IF(NTHINC.NE.1) GO TO 1773
      DFS=DTIME*ABS(DADOTG(KN)) 
      DFILMT(KN)=DFS/DRATLS     
      GO TO 1774      
 1773 CONTINUE
      ALPHA=TKOX/DCRUST(KN)     
      GAMMA=(DENOX*COX*DCRUST(KN)*TKF*(TFREEZ-TFRZBM))/(2.D0*DENS*DRATLS
     1*HABL)
      CON5=(TKF/ALPHA)*(TKF/ALPHA)+(4.D0*GAMMA)/ALPHA
      DFILMT(KN)=-0.5D0*(TKF/ALPHA)+0.5D0*SQRT(CON5)
      DFS=DFILMT(KN)*DRATLS     
 1774 CONTINUE
      IF(NTYPMT(IP,KN).GT.3) GO TO 3418   
      IF(NL(IP,KN).EQ.2) GO TO 3418      
      DO 3419 KTP=1,16
      DFILM(KTP,KN)=0.D0
      IF(KTP.LT.16) GO TO 3419  
      DFILM(KTP,KN)=DFILMT(KN)  
 3419 CONTINUE
      GO TO 3421      
 3418 CONTINUE
      DO 3420 KTP=1,16
      CALL INDEX(KTP,KJ)
      DENR=ROM(KJ)/ROMLIQ(KJ)   
      IF(KTP.EQ.16) DENR=RSLAGS/RSLAGL    
      DINCR=(DFS*HBMT(IP,KN,KTP))/ZLENGT  
      IF(NL(IP,KN).EQ.2) DINCR=XWTSS(KTP)*DFS      
      DFILM(KTP,KN)=DENR*DINCR  
      IF(NL(IP,KN).EQ.2) GO TO 3420      
      HBMT(IP,KN,KTP)=HBMT(IP,KN,KTP)-DINCR
      HBMT(IP,KN,KTP)=DMAX1(HBMT(IP,KN,KTP),0.D0)   
 3420 CONTINUE
 3421 CONTINUE
      ZABLAT(KN)=ZABLAT(KN)-DFS 
      ZABOLD(KN)=ZABLAT(KN)
      ELEVAT(KN)=ELEVAT(KN)+DCRUST(KN)-DFS+DFILMT(KN)
      IF(NTYPMT(1,KN).LE.3) DABCON(KN)=DABOLD(KN)+DFS
      ELOLD(KN)=ELEVAT(KN)
      DABOLD(KN)=DABCON(KN)
      GO TO 7750      
C.....CHECK THERMAL CONDITIONS AT MELT (OR CRUST)/BASEMAT INTERFACE TO  
C.....DETERMINE IF ABLATION CAN ENSUE AT THIS INTERFACE.      
 7751 CONTINUE
      IF(NPASS(KN).EQ.1.AND.NBCINT(KN).EQ.2) GO TO 7750      
      IF(NPASS(KN).EQ.1.AND.NBCINT(KN).EQ.3) GO TO 7750
      IF(NBCINT(KN).GE.4) GO TO 7750      
      IF(NTYPMT(1,KN).GT.3) GO TO 7752    
      TFRZCN=TSCS(NL(1,KN))     
      IF(TEMP(1,KN).LE.TFRZCN) GO TO 7750 
      ESOL=ESCS(NL(1,KN))
      HABL=ESCL(NL(1,KN))-ESCS(NL(1,KN))  
      DTDES=1.D0/CCSS(NL(1,KN)) 
      DRAT=RMASSS/RSLAGL
      IF(NL(1,KN).EQ.2) DRAT=ROSTLS/ROSTLL
      GO TO 7753      
 7752 CONTINUE
      IF(NTYPMT(1,KN).EQ.5) GO TO 7754    
      TFRZL=TFML      
      TFRZCN=TFMS     
      GO TO 7755      
 7754 CONTINUE
      TFRZL=TFOL      
      TFRZCN=TFOS     
 7755 CONTINUE
      IF(TEMP(1,KN).LE.TFRZCN) GO TO 7750 
      CALL ASINED(1,KN,HSEND,ROMLIQ,WPCC,WPM,WPA,WPCS,ROC,RSLAGL,0)     
      CALL ETF(EFL,HSEND,TFRZL,DTDE)      
      CALL ASINED(1,KN,HSEND,ROM,WPCC,WPM,WPA,WPCS,ROC,RSLAGS,0)
      CALL ETF(ESOL,HSEND,TFRZCN,DTDES)   
      HABL=EFL-ESOL   
      CALL MASFRD(ROM,1,KN,RSLAGS,ROC,DENS)
      CALL MASFRD(ROMLIQ,1,KN,RSLAGL,ROC,DENL)      
      DRAT=DENS/DENL  
 7753 CONTINUE
C.....ALGORITHM REACHES HERE, BASEMAT ABLATION STARTS.  IF INITIALLY    
C.....NBC=1 SET NBC=4 OR 5 DEPENDING ON THE SETTING OF NABLFM.  IF      
C.....INITIALLY NBC=2 OR 3 SET NBC=6 OR 8 DEPENDING UPON THE SETTING OF 
C.....NDRNFM.
      IF(NBCINT(KN).GT.1) GO TO 6744      
      IF(NABLFM.EQ.1) GO TO 7756
      NBCINT(KN)=4    
      GO TO 7757      
 7756 CONTINUE
      IF(HTOT(KN).LE.0.D0) GO TO 7750     
      NBCINT(KN)=5    
      GO TO 7757      
 6744 CONTINUE
      IF(NDRNFM.EQ.1) GO TO 6746
      NBCINT(KN)=8    
      IF(NTHINC.EQ.1) NBCINT(KN)=9
      GO TO 7757      
 6746 CONTINUE
      NBCINT(KN)=6    
      IF(NTHINC.EQ.1) NBCINT(KN)=7
 7757 CONTINUE
      DFS=(XBTW(1,KN)*(ENTHP(1,KN)-ESOL))/HABL      
      DFSMAX=XBTW(1,KN)*(1.D0-XFCABL)     
      DFS=DMIN1(DFS,DFSMAX)     
      DFM=DFS*DRAT    
      IF(NBCINT(KN).EQ.4) GO TO 7758      
      IF(NBCINT(KN).EQ.6.OR.NBCINT(KN).EQ.7) GO TO 7758
      GO TO 7759      
 7758 CONTINUE
C.....ADD INITIALLY ABLATED BASEMAT MATERIAL TO BULK MELT.    
      IF(HTOT(KN).GT.0.D0) GO TO 8818     
      TBULK(KN)=TFRZCN
      ENBLK(KN)=ESOL  
 8818 CONTINUE
      HTOT(KN)=HTOT(KN)+DFM     
      IF(NTYPMT(1,KN).GT.3) GO TO 7760    
      IF(NL(1,KN).EQ.2) GO TO 7760
      HOXIDE(KN)=HOXIDE(KN)+DFM 
      HITE(16,KN)=HITE(16,KN)+DFM
      GO TO 7762      
 7760 CONTINUE
      HMETAL(KN)=0.D0 
      HOXIDE(KN)=0.D0 
      DO 7763 KTP=1,16
      CALL INDEX(KTP,KJ)
      DENR=ROM(KJ)/ROMLIQ(KJ)   
      IF(KTP.EQ.16) DENR=RSLAGS/RSLAGL    
      DINCR=(DFS*HBMT(1,KN,KTP))/XBTW(1,KN)
      IF(NL(1,KN).EQ.2) DINCR=XWTSS(KTP)*DFS
      HITE(KTP,KN)=HITE(KTP,KN)+DINCR*DENR
      IF(KTP.LE.7) HMETAL(KN)=HMETAL(KN)+HITE(KTP,KN)
      IF(KTP.GE.8) HOXIDE(KN)=HOXIDE(KN)+HITE(KTP,KN)
      IF(NL(1,KN).EQ.2) GO TO 7763
      HBMT(1,KN,KTP)=HBMT(1,KN,KTP)-DINCR 
      HBMT(1,KN,KTP)=DMAX1(HBMT(1,KN,KTP),0.D0)     
 7763 CONTINUE
      GO TO 7762      
 7759 CONTINUE
C.....NBC=5 OR 8 CASES TREATED HERE.  PLACE INITALLY ABLATED MATERIAL   
C.....INTO MOLTEN FILM.
      DFILMT(KN)=DFM  
      IF(NTYPMT(1,KN).GT.3) GO TO 7764    
      IF(NL(1,KN).EQ.2) GO TO 7764
      DO 7765 KTP=1,16
      DFILM(KTP,KN)=0.D0
      IF(KTP.LT.16) GO TO 7765  
      DFILM(KTP,KN)=DFILMT(KN)  
 7765 CONTINUE
      GO TO 7762      
 7764 CONTINUE
      DO 7766 KTP=1,16
      CALL INDEX(KTP,KJ)
      DENR=ROM(KJ)/ROMLIQ(KJ)   
      IF(KTP.EQ.16) DENR=RSLAGS/RSLAGL    
      DINCR=(DFS*HBMT(1,KN,KTP))/XBTW(1,KN)
      IF(NL(1,KN).EQ.2) DINCR=XWTSS(KTP)*DFS
      DFILM(KTP,KN)=DENR*DINCR  
      IF(NL(1,KN).EQ.2) GO TO 7766
      HBMT(1,KN,KTP)=HBMT(1,KN,KTP)-DINCR 
 7766 CONTINUE
 7762 CONTINUE
C.....FOR ABLATION CASES, SET RELEVANT INFORMATION FOR START OF ABLATION
C.....CALCULATION.    
      TEMP(1,KN)=TFRZCN
      ENTHP(1,KN)=ESOL
      DTDEB(1,KN)=DTDES
      IF(NTYPMT(1,KN).LE.3.AND.NL(1,KN).EQ.1) EMAX(1,KN)=ENTHP(1,KN)    
      ZABLAT(KN)=XDIST(1,KN)-DFS
      ZABOLD(KN)=ZABLAT(KN)
      XDISTO(KN)=XDIST(1,KN)    
      XBTWO(KN)=XBTW(1,KN)      
      IF(NBCINT(KN).EQ.4.OR.NBCINT(KN).EQ.6.OR.NBCINT(KN).EQ.7) ELEVAT(K
     1N)=ELEVAT(KN)-DFS
      IF(NBCINT(KN).EQ.5.OR.NBCINT(KN).EQ.8.OR.NBCINT(KN).EQ.9) ELEVAT(K
     1N)=ELEVAT(KN)-DFS+DFM     
      IF(NTYPMT(1,KN).LE.3) DABCON(KN)=DABOLD(KN)+DFS
      ELOLD(KN)=ELEVAT(KN)
      DABOLD(KN)=DABCON(KN)
 7750 CONTINUE
C.....IF NEEDED, RECOMPUTE BULK MELT GLOBAL INVENTORIES WHEN CRUSTING OR 
C.....ABLATION WITH DRAINAGE TO MELT HAS BEEN INITIATED AT THIS TIMESTEP.
      IF(NCHNG.EQ.0) GO TO 88
      TOTMET=0.D0     
      TOTOX=0.D0      
      TOTVOL=0.D0     
      DO 3447 KT=1,16 
      VOLCN(KT)=0.D0  
      DO 3449 KN=1,NUMNOD
      VOLCN(KT)=VOLCN(KT)+HITE(KT,KN)*AREA(KN)*FLOAT(NDOOR)   
 3449 CONTINUE
      TOTVOL=TOTVOL+VOLCN(KT)   
      IF(KT.LE.7) TOTMET=TOTMET+VOLCN(KT) 
      IF(KT.GE.8) TOTOX=TOTOX+VOLCN(KT)   
 3447 CONTINUE
      DO 7491 KN=1,NUMNOD
      HTOT(KN)=0.D0   
      DO 7492 KT=1,16 
      HTOT(KN)=HTOT(KN)+HITE(KT,KN)
 7492 CONTINUE
 7491 CONTINUE
   88 CONTINUE
      RETURN
      END