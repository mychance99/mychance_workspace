C...................................................................... 
C.....SUBROUTINE SHELL CALCULATES THE DRYWELL SHELL FORWARD ELIMINATION 
C.....COEFFICIENTS FOR THE SOLUTION OF SHELL ENTHALPIES & TEMPERATURES.  
C...................................................................... 
      SUBROUTINE SHELL(IFAILJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL(8) KL
      DIMENSION DENE(10,999),FKE(10,999),ALPE(10,999),TARA(999)   
      DIMENSION TARB(999),P(10),Q(10),A(10),B(10),C(10),D(10)  
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
      COMMON /WATVARS/FCRUST(999),VLW(999),VLWOLD(999),DWAT(999),
     1DWATOLD(999),EWAT(999),EWATOLD(999),TWAT(999),CORDDC(999),
     2DSRDC(999),DHDC(999),TSRDC,TINTDC,HWATB(999),TSURFW(999),XMS(999),
     3XMST,XINTS,ESRDTW(999),TSRDTW(999),DSRDTW(999),TSRDW,TINTDW,
     4TWATI,ELDCO(999),ELDCX(999,99),CRDCX(999,99),TIMINJ(999,99),
     5XDTINJ(999,99),TDTINJ(999,99),ELWATI,XMWAT,XMWATO,XBALW,
     6HDRY(999),EINTW,XMCT(199),XMCDT(16,999),TMPCDT(999),CMPCDT(999),
     7FDC(999),FCOV(999),FBED(999),FHEAT(999),XFACJ(999),
     8XMBEDJ(16,999),XMBEDJT(999),PBED,QBED(999),QWATER(999),
     9QSURFACE(999),QWATERT(999),HBED(999),POROSBED,XMCRDT(16,999),
     1XMCRT(999),XMBDINT(999),XMT(999)
      COMMON/PROPM/DTDEM(999),CPMLT(999),FKMLT(999),DENMLT(999),
     1SIGMLT(999),UMMLT(999),EMIMLT(999)      
      COMMON/DATSIMS/TSIS,TSIL,CSIS,CSIL,DHSIL,RSIS,RSIL,TKSIS,TKSIL,
     1ESIL,XCSI(16)
      DATA EMISS / 0.3D0 /    
C.....FIRST CHECK TO SEE IF THE CALCULATION IS TO BE PERORMED...
      IF(NSKIPE.EQ.0.OR.NADAB.EQ.1) GO TO 8017
      IF(ISHELE.EQ.1) GO TO 8017
      IF(NACSH.EQ.0) GO TO 8017 
C.....FIND HEAT TRANSFER COEFFICIENT TO SOLIDIFIED DEBRIS IF BULK
C.....FREEZING HAS OCCURED ADJACENT TO THE SHELL.   
      IF(NFRSCT.EQ.1) GO TO 1111
      HDEBRS=0.D0     
      GO TO 1112      
 1111 CONTINUE
      TEV=TIME-TFRZSH 
      IF(TEV.LT.DTIME) TEV=DTIME
      DENOM=PDEBR*CPDEBR
      AFDEBR=0.D0     
      IF(DENOM.GT.0.D0) AFDEBR=TKDEBR/DENOM
      DENOM=0.D0      
      IF(AFDEBR.GT.0.D0) DENOM=1.D0/SQRT(PI*TEV*AFDEBR)
      HDEBRS=TKDEBR*DENOM
 1112 CONTINUE
      IF(HTOT(NBMADJ).LE.0.D0) GO TO 9310 
      NS=NBMADJ
      ALPBLK=DTDEM(NS)
      CMS=CPMLT(NS)   
      PMS=DENMLT(NS)  
      UMS=UMMLT(NS)   
      TKMS=FKMLT(NS)  
      CALL SHELLM(TBULK(NS),VGJ(NS),VEL(NS),VOID(NS),CMS,PMS,TKMS,UMS,HT
     1OT(NS),HSHELL,XLENSH)     
      GO TO 9313      
 9310 CONTINUE
      HSHELL=0.D0     
      ALPBLK=0.D0     
 9313 CONTINUE
C.....FIND SHELL FORWARD ELIMINATION COEFFICIENTS.  FIRST FIND SHELL    
C.....PROPERTIES.     
      DO 9299 KV=1,NUMSHV
      DO 9300 KH=1,NUMSHH
      CALL DENSS(TENDP(KH,KV),DENE(KH,KV))
      CALL CONDSS(TENDP(KH,KV),FKE(KH,KV))
      CALL TESS(ENENDP(KH,KV),TPS,ALPE(KH,KV))      
 9300 CONTINUE
 9299 CONTINUE
      ANFAR=(PI*ANGSHL)/180.D0  
      ANSIN=COS(ANFAR)
      VFMAX=VOID(NBMADJ)
      IF(VFMAX.GE.0.3) VFMAX=0.3
      THM=HTOT(NBMADJ)/(ANSIN*(1.D0-VFMAX))+ELEVAT(NBMADJ)/ANSIN 
      IF(THM.LT.HNODT) GO TO 8830 
      WRITE(6,7671) TIME,HNODT,HTOT(NBMADJ),THM
      WRITE(2,7671) TIME,HNODT,HTOT(NBMADJ),THM
 7671 FORMAT(//,'****MELT DEPTH ADJACENT TO SHELL HAS EXCEEDED THE NODAL
     1IZED LENGTH; INCREASE HNODT AND RESTART',/,' TIME=',G12.5,/,' NODA
     2LIZED LENGTH=',G12.5,/,' MELT DEPTH=',G12.5,/,' LENGTH COVERED BY 
     3VOIDED MELT=',G12.5)
      IFAILJ=1
      GO TO 110
 8830 CONTINUE 
      HENP=-DXVERT/2.D0
      IF(NWAT.GT.0) CALL CONWAT(TSAT,KL,PL,CL,UL,HLV,EWL,EWV,CWV,SIGMA,
     1PDRYWL)
      TSUB=0.D0
      IF(NWAT.EQ.2) TSUB=DMAX1(TSAT-TWAT(NBMADJ),0.D0)
      DO 9301 KV=1,NUMSHV
C.....ASSIGN CRUST PROPERTIES.  
      TKC=FKOXE(KV)   
      SIGC=SIGOXE(KV) 
      IF(KV.EQ.1) GO TO 8924    
      HENP=HENP+DXVERT
      HENT=HENP+DXVERT/2.D0     
      HENB=HENP-DXVERT/2.D0     
 8924 CONTINUE
      XVERT=DXSNK     
      IF(KV.GT.1) XVERT=DXVERT  
      DO 9319 KH=1,NUMSHH
C.....CALCULATE CONDUCTANCES.   
      CHM=0.D0
      IF(KH.NE.1) CHM=(FKE(KH-1,KV)*FKE(KH,KV))/(XBTE(KH-1,KV)*FKE(KH,KV
     1)+(XDSTE(KH-1,KV)-XBTE(KH-1,KV))*FKE(KH-1,KV))
      CHP=0.D0
      IF(KH.NE.NUMSHH) CHP=(FKE(KH,KV)*FKE(KH+1,KV))/(XBTE(KH,KV)*FKE(KH
     1+1,KV)+(XDSTE(KH,KV)-XBTE(KH,KV))*FKE(KH,KV)) 
      CVM=0.D0
      IF(KV.EQ.2) CVM=(2.D0*FKE(KH,KV-1)*FKE(KH,KV))/(DXSNK*FKE(KH,KV)+D
     1XVERT*FKE(KH,KV-1))
      IF(KV.GT.2) CVM=(2.D0*FKE(KH,KV-1)*FKE(KH,KV))/(DXVERT*(FKE(KH,KV)
     1+FKE(KH,KV-1))) 
      CVP=0.D0
      IF(KV.EQ.1) CVP=(2.D0*FKE(KH,KV)*FKE(KH,KV+1))/(DXSNK*FKE(KH,KV+1)
     1+DXVERT*FKE(KH,KV))
      IF(KV.GT.1) CVP=(2.D0*FKE(KH,KV)*FKE(KH,KV+1))/(DXVERT*(FKE(KH,KV)
     1+FKE(KH,KV+1))) 
      IF(KH.GT.1) GO TO 9302    
      IF(KV.GT.1) GO TO 9303    
      HXLA(KV)=0.D0   
      HXLB(KV)=0.D0   
      TARA(KV)=0.D0   
      TARB(KV)=0.D0   
      TSFEB(KV)=0.D0  
      GO TO 9304      
 9303 CONTINUE
      XFACM=(THM-HENB)/(HENT-HENB)
      XFACM=DMAX1(XFACM,0.D0)   
      IF(XFACM.GT.1.D0) XFACM=1.D0
      XFACA=1.D0-XFACM
      IF(ELEVAT(NBMADJ)/ANSIN.LT.HENB) GO TO 1113   
      TARA(KV)=TDEBRS 
      HXLA(KV)=HDEBRS*DXVERT    
      HXLB(KV)=0.D0   
      TARB(KV)=0.D0   
      TSFEB(KV)=TENDP(KH,KV)    
      GO TO 9304      
 1113 CONTINUE
      TSFEB(KV)=TENDP(KH,KV)    
      IF(NBSHL(KV).EQ.1) GO TO 9400
      IF(NBSHL(KV).EQ.3) GO TO 9400
      TSFEB(KV)=TEFZX(KV)
 9400 CONTINUE
      IF(XFACA.GT.0.D0) GO TO 9401
      HXLA(KV)=0.D0   
      TARA(KV)=0.D0   
      GO TO 9402      
 9401 CONTINUE
      IF(NWAT.GT.0) GO TO 9403  
      TARA(KV)=TBOUND 
      RADCST=STEF/(1.D0/EMISCN+1.D0/EMISS-1.D0)     
      HTC=RADCST*(TSFEB(KV)**2+TARA(KV)**2)*(TSFEB(KV)+TARA(KV))
      IF(NBSHL(KV).EQ.3) HTC=(HTC*TKC)/(DCRS(KV)*HTC+
     1TKC)  
      HXLA(KV)=XFACA*DXVERT*HTC 
      GO TO 9402      
 9403 CONTINUE
      TARA(KV)=TSAT   
      CALL SHELLW(TSFEB(KV),TSUB,PDRYWL,HTC,ANGSHL)      
      IF(NBSHL(KV).EQ.3) HTC=(HTC*TKC)/(DCRS(KV)*HTC+
     1TKC)  
      HXLA(KV)=XFACA*DXVERT*HTC 
      IF(TSFEB(KV).LT.TSAT+1.0) GO TO 666
  666 CONTINUE
 9402 CONTINUE
      IF(XFACM.GT.0.D0) GO TO 9405
      TARB(KV)=0.D0   
      HXLB(KV)=0.D0   
      GO TO 9304      
 9405 CONTINUE
      TARB(KV)=TBULK(NBMADJ)    
      HXLB(KV)=XFACM*DXVERT*HSHELL
 9304 CONTINUE
      A1=(DENE(KH,KV)*XVERT*XBTE(KH,KV))/DTIME      
      A2=XBTE(KH,KV)*ALPE(KH,KV)*(CVP+CVM)
      IF(NBSHL(KV).NE.2) GO TO 4353  
      IF(DCRS(KV).GT.0.D0) GO TO 9353  
 4353 CONTINUE
      IF(NBSHL(KV).NE.3) GO TO 7322 
      IF(XFACM.GT.0.D0.AND.DCRS(KV).GT.0.D0) GO TO 9351
 7322 CONTINUE
      A3=CHP*XVERT*ALPE(KH,KV)+(HXLA(KV)+HXLB(KV))*ALPE(KH,KV)
      GO TO 9352      
 9351 CONTINUE
      A3=XVERT*(TKC/DCRS(KV)+CHP)*ALPE(KH,KV)
      GO TO 9352      
 9353 CONTINUE
      HTCC=HXLB(KV)/XVERT
      SG1=(HTCC*(TBULK(NBMADJ)-TSFEB(KV)))/SIGC     
      DTEV=TSFEB(KV)-TENDP(KH,KV)
      DTEV=DMAX1(DTEV,0.D0)     
      SG2=4.D0*(DCRSLD(KV)**2+(TKC*DTEV)/SIGC)      
      GAM=SQRT(SG1**2+SG2)      
      WJP1=-TKC/(SIGC*GAM)      
      SB1=((TBULK(NBMADJ)-TSFEB(KV))*HTCC**2)/(GAM*SIGC**2)   
      WBLK=0.5D0*(SB1-HTCC/SIGC)
      A3=XVERT*(CHP+(TKC*(TSFEB(KV)-TENDP(KH,KV))*WJP1)/DCRS(KV)**2+TKC/
     1DCRS(KV))*ALPE(KH,KV)     
 9352 CONTINUE
      A(KH)=A1+A2+A3  
      B(KH)=XVERT*CHP*ALPE(KH+1,KV)
      C(KH)=0.D0      
      D1=(DENE(KH,KV)*XVERT*XBTE(KH,KV)*(ENENDP(KH,KV)-ENOLDP(KH,KV)))/D
     1TIME  
      D2=0.D0
      IF(KH.LT.NUMSHH) D2=XVERT*CHP*(TENDP(KH,KV)-TENDP(KH+1,KV)) 
      D3=0.D0
      IF(KV.GT.1) D3=XBTE(KH,KV)*CVM*(TENDP(KH,KV-1)-TENDP(KH,KV)) 
      D4=0.
      IF(KV.LT.NUMSHV) D4=XBTE(KH,KV)*CVP*(TENDP(KH,KV)-TENDP(KH,KV+1))
      IF(NBSHL(KV).NE.2) GO TO 6211
      IF(DCRS(KV).GT.0.D0) GO TO 9355
 6211 CONTINUE
      IF(NBSHL(KV).NE.3) GO TO 7211
      IF(XFACM.GT.0.D0.AND.DCRS(KV).GT.0.D0) GO TO 9355 
 7211 CONTINUE
      D5=HXLA(KV)*(TARA(KV)-TENDP(KH,KV))+HXLB(KV)*(TARB(KV)-TENDP(KH,KV
     1))    
      OMEG=HXLB(KV)*ALPBLK      
      GO TO 9357      
 9355 CONTINUE
      D5=(XVERT*TKC*(TSFEB(KV)-TENDP(KH,KV)))/DCRS(KV)
      OMEG=0.D0
      IF(NBSHL(KV).EQ.2) OMEG=-(D5*WBLK*ALPBLK)/DCRS(KV)      
 9357 CONTINUE
      D(KH)=-D1-D2+D3-D4+D5     
      GO TO 9319      
 9302 CONTINUE
      XHOR=XDSTE(KH-1,KV)-XBTE(KH-1,KV)   
      IF(KH.NE.NUMSHH) XHOR=XBTE(KH,KV)+XHOR
      A1=(DENE(KH,KV)*XVERT*XHOR)/DTIME   
      A2=XVERT*(CHM+CHP)*ALPE(KH,KV)      
      A3=XHOR*(CVP+CVM)*ALPE(KH,KV)
      A(KH)=A1+A2+A3  
      B(KH)=0.
      IF(KH.LT.NUMSHH) B(KH)=XVERT*CHP*ALPE(KH+1,KV)
      C(KH)=XVERT*CHM*ALPE(KH-1,KV)
      D1=(DENE(KH,KV)*XVERT*XHOR*(ENENDP(KH,KV)-ENOLDP(KH,KV)))/DTIME   
      D2=XVERT*CHM*(TENDP(KH-1,KV)-TENDP(KH,KV))    
      D3=0.
      IF(KH.LT.NUMSHH) D3=XVERT*CHP*(TENDP(KH,KV)-TENDP(KH+1,KV))     
      D4=0.
      IF(KV.GT.1) D4=XHOR*CVM*(TENDP(KH,KV-1)-TENDP(KH,KV))
      D5=0.
      IF(KV.LT.NUMSHV) D5=XHOR*CVP*(TENDP(KH,KV)-TENDP(KH,KV+1))     
      D(KH)=-D1+D2-D3+D4-D5     
 9319 CONTINUE
C.....CONSTRUCT SHELL FORWARD ELIMINATION COEFFICIENTS GIVEN TDMA MATRIX
C.....FOR THE KH'TH HORIZONTAL NODAL VECTOR.
      P(1)=B(1)/A(1)  
      Q(1)=D(1)/A(1)  
      DO 9325 KD=2,NUMSHH
      P(KD)=B(KD)/(A(KD)-C(KD)*P(KD-1))   
      Q(KD)=(D(KD)+C(KD)*Q(KD-1))/(A(KD)-C(KD)*P(KD-1))
 9325 CONTINUE
      THETE0(NUMSHH,KV)=Q(NUMSHH)
      DO 9326 KD=1,NUMSHH-1     
      IARG=NUMSHH-KD  
      THETE0(IARG,KV)=P(IARG)*THETE0(IARG+1,KV)+Q(IARG)
 9326 CONTINUE
      IF(ALPBLK.GT.0.D0.AND.HXLB(KV).GT.0.D0) GO TO 9327      
      DO 9328 KD=1,NUMSHH
      THETE1(KD,KV)=0.D0
 9328 CONTINUE
      GO TO 9301      
 9327 CONTINUE
      Q(1)=OMEG/A(1)  
      DO 9329 KD=2,NUMSHH
      Q(KD)=(C(KD)*Q(KD-1))/(A(KD)-C(KD)*P(KD-1))   
 9329 CONTINUE
      THETE1(NUMSHH,KV)=Q(NUMSHH)
      DO 9330 KD=1,NUMSHH-1     
      IARG=NUMSHH-KD  
      THETE1(IARG,KV)=P(IARG)*THETE1(IARG+1,KV)+Q(IARG)
 9330 CONTINUE
 9301 CONTINUE
C.....DETERMINE INTEGRATED HEAT LOAD INTO SHELL FOR SOLUTION OF BULK    
C.....CONSERVATION OF ENERGY EQUATION.    
      QSHELL=0.D0     
      QSHELE=0.D0     
      DO 9398 KD=1,NUMSHV
      QSHELL=QSHELL+HXLB(KD)*(TBULK(NBMADJ)-TSFEB(KD)-ALPE(1,KD)*THETE0(
     11,KD))
      QSHELE=QSHELE+HXLB(KD)*(ALPBLK-ALPE(1,KD)*THETE1(1,KD)) 
 9398 CONTINUE
 8017 CONTINUE
  110 CONTINUE
      RETURN
      END