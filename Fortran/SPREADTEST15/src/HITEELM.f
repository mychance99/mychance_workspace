C.......................................................................
C.....SUBROUTINE HITEELM CALCULATES BOTH THE HEIGHT INCREMENTS AND 
C.....HEIGHTS FOR THE COM EQUATION DEPENDING UPON THE SETTING OF NWRAP
C.....WHICH IS SET IN THE MAIN PROGRAM.
C.......................................................................
      SUBROUTINE HITEELM(NWRAP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
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
      COMMON/BASPR/ HBMT(999,999,16),HBMOLD(999,999,16) 
C.....FIRST ZERO OUT HEIGHT INCEMENTS FROM PREVIOUS TIMESTEP.  
      DO 4064 KN=1,NUMNOD
      DO 4063 KTP=1,16
      DELHIT(KTP,KN)=0.D0
 4063 CONTINUE
 4064 CONTINUE
C.....NOW SOLVE THE TDMA FOR HEIGHT INCREMENTS OR HEIGHT, DEPENDING
C.....UPON NWRAP.
      DO 5063 KTP=1,16
C.....CHECK IF KTP-TH MELT CONSTITUENT HAS BEEN EJECTED FROM THE RPV.   
      IF(NTRACK(KTP).EQ.0) GO TO 5087     
      CALL INDEX(KTP,KJIM)      
      DO 5064 KND=1,NUMNOD      
C.....SET UP THE TDMA PROBLEM FOR THE KTP-TH MELT CONSTITUENT.  FIRST   
C.....DETERMINE THE OXIDATION CONTANTS FOR ZIRCONIUM, IRON, AND CHROMIUM
C.....LAYERS.
      ZCOMP=XDISTO(KND)-ZABLAT(KND)
      IF(ZCOMP.GE.XBTWO(KND)) GO TO 5081  
      NBM=1 
      XLENGT=ZABLAT(KND)-(XDISTO(KND)-XBTWO(KND))   
      GO TO 5082      
 5081 CONTINUE
      NBM=2 
      XLENGT=ZABLAT(KND)+XBTW(2,KND)      
 5082 CONTINUE
      IF(HTMP(KND).LE.0.D0) GO TO 5069    
      CALL INDEX(KTP,KJIM)      
      IF(KTP.LE.3) GO TO 5065   
      IF(KTP.EQ.8) GO TO 5066   
      IF(KTP.EQ.9) GO TO 5067   
      IF(KTP.EQ.12) GO TO 5068  
 5069 CONTINUE
      GH2ORE=0.D0     
      GCO2RE=0.D0     
      GO TO 5070      
 5065 CONTINUE
      IF(HITE(KTP,KND).LE.0.D0) GO TO 5069
      GO TO(5071,5072,5073),KTP 
 5071 CONTINUE
      GCON=FMMOL(KJIM)/(ROMLIQ(KJIM)*XZRMRE)
      GO TO 5074      
 5072 CONTINUE
      IF(HITE(1,KND).GT.0.D0) GO TO 5069  
      IF(HITE(3,KND).GT.0.D0) GO TO 5069  
      GCON=FMMOL(KJIM)/(3.D0*ROMLIQ(KJIM)*XFEMRE)   
      GO TO 5074      
 5073 CONTINUE
      IF(HITE(1,KND).GT.0.D0) GO TO 5069  
      GCON=FMMOL(KJIM)/(ROMLIQ(KJIM)*XCRMRE)
      GO TO 5074      
 5066 CONTINUE
      IF(HITE(1,KND).LE.0.D0) GO TO 5069  
      GCON=-FMMOL(KJIM)/(ROMLIQ(KJIM)*XZRORE)
      GO TO 5074      
 5067 CONTINUE
      IF(HITE(1,KND).GT.0.D0) GO TO 5069  
      IF(HITE(3,KND).GT.0.D0) GO TO 5069  
      IF(HITE(2,KND).LE.0.D0) GO TO 5069  
      GCON=-FMMOL(KJIM)/(ROMLIQ(KJIM)*XFEORE)
      GO TO 5074      
 5068 CONTINUE
      IF(HITE(1,KND).GT.0.D0) GO TO 5069  
      IF(HITE(3,KND).LE.0.D0) GO TO 5069  
      GCON=-FMMOL(KJIM)/(ROMLIQ(KJIM)*XCRORE)
 5074 CONTINUE
      TSEND=TEMP(1,KND)
      IF(TSEND.GT.TCS) TSEND=TCS
      DBSEND=DMIN1(HTOT(KND),DBUBOX)      
      DBSEND=DMAX1(HMINC,DBSEND)
      CALL OXID(PDRYWL,TSEND,TBULK(KND),DBSEND,HTOT(KND),UTRISE,OXH2O,OX
     1CO2)  
      GH2ORE=(AREA(KND)*GCON*OXH2O)/XMH2O 
      GCO2RE=(AREA(KND)*GCON*OXCO2)/XMCO2 
      IF(KTP.EQ.2.OR.KTP.EQ.9) GCO2RE=0.D0
 5070 CONTINUE
C.....CALCULATE CRUSTING AND/OR ABLATION CONSTANTS FOR THE VARIOUS      
C.....TYPES OF MELT/CONCRETE INTERFACIAL BOUNDARY CONDITIONS. 
      IF(HTMP(KND).LE.0.D0) GO TO 8311    
      IF(NBCINT(KND).EQ.2.AND.NCRRT(KND).EQ.1) GO TO 5075     
      IF(NBCINT(KND).EQ.6.AND.NCRRT(KND).EQ.1) GO TO 5075     
      IF(NBCINT(KND).EQ.8.AND.NCRRT(KND).EQ.1) GO TO 5075     
 8311 CONTINUE
      THETCR=0.D0     
      GO TO 5076      
 5075 CONTINUE
C.....CRUSTING CONSTANTS ARE DETERMINED BELOW.      
      IF(NCRSTM(KND).EQ.1.AND.KTP.LE.7) GO TO 8311  
      IF(NCRSTM(KND).EQ.2.AND.KTP.GE.8) GO TO 8311  
      IF(HITE(KTP,KND).LE.0.D0) GO TO 8311
      RDRAT=ROMLIQ(KJIM)/ROM(KJIM)
      IF(KTP.EQ.16) RDRAT=RSLAGL/RSLAGS   
      THETCR=(RDRAT*AREA(KND)*HITE(KTP,KND))/HTOT(KND)
 5076 CONTINUE
      IF(HTMP(KND).LE.0.D0) GO TO 5079    
      IF(NBCINT(KND).EQ.4) GO TO 5077     
      IF(NBCINT(KND).EQ.6) GO TO 5077     
      IF(NBCINT(KND).EQ.7) GO TO 5077     
 5079 CONTINUE
      THETAB=0.D0     
      GO TO 5078      
 5077 CONTINUE
C.....ABLATION CONSTANTS ARE DETERMINED BELOW.  NOTE THAT THE PROPER    
C.....BASEMAT PROPERTY INDICES HAVE TO BE SET DEPENDING ON PENETRATION  
C.....DEPTH OF THE ABLATION FRONT.
      IF(NTYPMT(NBM,KND).GT.3) GO TO 6078 
      IF(NL(NBM,KND).EQ.2) GO TO 6080    
      IF(KTP.LT.16) GO TO 5079  
C.....CONCRETE ABLATION CONSTANT ASSIGNED BELOW.    
      THETAB=(RMASSS*AREA(KND))/RSLAGL    
      GO TO 5078      
C.....FOR BASEMATS WHICH ARE COMPOSED NOT OF CONCRETE, THE FOLLOWING    
C.....DETERMINES THE ABLATION CONTANTS.   
 6080 CONTINUE
      THETAB=(ROM(KJIM)*AREA(KND)*XWTSS(KTP))/ROMLIQ(KJIM)    
      GO TO 5078      
 6078 CONTINUE
      IF(HBMT(NBM,KND,KTP).LE.0.D0) GO TO 5079      
      RDRAT=ROM(KJIM)/ROMLIQ(KJIM)
      IF(KTP.EQ.16) RDRAT=RSLAGS/RSLAGL   
      THETAB=(RDRAT*AREA(KND)*HBMT(NBM,KND,KTP))/XLENGT
 5078 CONTINUE
C.....ASSIGN VALUES FOR THE EVALUATION OF THE TDMA ELEMENTS.  
      VJM12=VEL(KND)  
      VJP12=VEL(KND+1)
      IF(KND.EQ.1) GO TO 5083   
      IF(KND.EQ.NUMNOD) GO TO 5084
      HJM1=HITE(KTP,KND-1)      
      HJ=HITE(KTP,KND)
      HJP1=HITE(KTP,KND+1)      
      GO TO 5085      
 5083 CONTINUE
      HJM1=0.D0
      HJ=HITE(KTP,KND)
      HJP1=HITE(KTP,KND+1)      
      GO TO 5085      
 5084 CONTINUE
      HJM1=HITE(KTP,KND-1)      
      HJ=HITE(KTP,KND)
      HJP1=0.D0
 5085 CONTINUE
C.....EVALUATE THE TDMA COEFFICIENTS.     
      IF(NWRAP.EQ.1) GO TO 5034 
      A1=0.D0
      GO TO 5035      
 5034 CONTINUE
      ABLK1=GH2ORE*(OMEG1(1,KND)+OMEG1(2,KND))+GCO2RE*(OMEG1(3,KND)+    
     1OMEG1(4,KND))+THETAB*TAO1(KND)+THETCR*SIG1(KND)
      A1=ABLK1*DBLKEN(KND)      
 5035 CONTINUE
      A2=(AREA(KND)*(HJ-HITOLD(KTP,KND)))/DTIME     
      A3=ARC(KND)*(HJ*DMAX1(-VJM12,0.D0)-HJM1*DMAX1(VJM12,0.D0))
      A4=ARC(KND+1)*(HJ*DMAX1(VJP12,0.D0)-HJP1*DMAX1(-VJP12,0.D0))      
      A5=GCO2RE*(SMFLX(3,KND)+SMFLX(4,KND))+GH2ORE*(SMFLX(1,KND)+SMFLX(2
     1,KND))+GH2ORE*(OMEG1(1,KND)+OMEG1(2,KND))+GCO2RE*(OMEG1(4,KND)    
     1+OMEG1(3,KND))  
      A6=THETAB*(DADOT(KND)+TAO0(KND))+THETCR*(DCDOT(KND)+SIG0(KND))    
      CALL INDEX(KTP,KJIM)      
      RODEN=ROMLIQ(KJIM)
      IF(KTP.EQ.16) RODEN=RSLAGL*ROC      
      A7=(AREA(KND)*SOURCE(KTP,KND))/RODEN
      D(KND)=-(A1+A2+A3+A4+A5+A6-A7)      
      C(KND)=ARC(KND)*DMAX1(VJM12,0.D0)   
      B(KND)=ARC(KND+1)*DMAX1(-VJP12,0.D0)
      A(KND)=AREA(KND)/DTIME+ARC(KND)*DMAX1(-VJM12,0.D0)+ARC(KND+1)*    
     1DMAX1(VJP12,0.D0)
 5064 CONTINUE
C.....EXECUTE TDMA ALGORITHM TO FIND THE INCREMENTAL HEIGHT INCREMENTS  
C.....FOR THE KTP-TH MELT CONSTITUENT.    
      P(1)=B(1)/A(1)  
      Q(1)=D(1)/A(1)  
      DO 5086 KTOP=2,NUMNOD     
      P(KTOP)=B(KTOP)/(A(KTOP)-C(KTOP)*P(KTOP-1))   
      Q(KTOP)=(D(KTOP)+C(KTOP)*Q(KTOP-1))/(A(KTOP)-C(KTOP)*P(KTOP-1))   
 5086 CONTINUE
C.....PERFORM BACKWARDS SUBSTITUTION TO FIND MELT CONSTITUENT HEIGHT    
C.....INCREMENTS.     
      DELHIT(KTP,NUMNOD)=Q(NUMNOD)
      DO 6159 KLP=1,NUMNOD-1    
      IARG=NUMNOD-KLP 
      DELHIT(KTP,IARG)=P(IARG)*DELHIT(KTP,IARG+1)+Q(IARG)     
 6159 CONTINUE
      GO TO 5063      
 5087 CONTINUE
      DO 5088 KNDP=1,NUMNOD     
      DELHIT(KTP,KNDP)=0.D0     
      HITE(KTP,KNDP)=0.D0
 5088 CONTINUE
 5063 CONTINUE
C.....GIVEN THE CONSTITUENT HEIGHT INCREMENTS, FIND THE TOTAL MELT      
C.....HEIGHT INCREMENT.
      DO 5022 KNDDR=1,NUMNOD    
      HTOT(KNDDR)=0.D0
      DELHTS(KNDDR)=0.D0
      DO 5023 KMT=1,16
      DELADD=DELHIT(KMT,KNDDR)
      IF(HITE(KMT,KNDDR)+DELADD.GE.0.D0) GO TO 7981
      DELADD=-HITE(KMT,KNDDR)
 7981 CONTINUE
      HCOMPI=HITE(KMT,KNDDR)+DELADD
      HITE(KMT,KNDDR)=DMAX1(HCOMPI,0.D0)  
      HTOT(KNDDR)=HTOT(KNDDR)+HITE(KMT,KNDDR)
      DELHTS(KNDDR)=DELHTS(KNDDR)+DELADD
 5023 CONTINUE
 5022 CONTINUE
      RETURN
      END