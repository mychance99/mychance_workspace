C..................................................................... 
C.....SUBROUTINE READER READS THE INPUT FILE FROM SPREAD.DAT AND DOES
C.....SOME MINOR INITIALIZATIONS IN ORDER TO START THE CALCULATION.
C..................................................................... 
      SUBROUTINE READER
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /LABEL/ TITLE(72)
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
      COMMON/DATSIMS/TSIS,TSIL,CSIS,CSIL,DHSIL,RSIS,RSIL,TKSIS,TKSIL,
     1ESIL,XCSI(16)
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
      COMMON /WATINT/NDOWNC,NDC(999),NPTDC(999),NINJ,NINJP(999),
     1NPTINJ(999)
C.....WRITE HEADER TO OUTPUT FILE FOR ECHO OF INPUT DATA.
      WRITE(6,781) 
  781 FORMAT(1X,'ECHO-WRITE OF INPUT DECK FOR THIS CALCULATION',//) 
C.....FIRST READ TEXT LABEL FOR INPUT FILE AND WRITE THAT BACK OUT TO 
C.....THE OUTPUT FILE.
      READ(5,6) TITLE
      WRITE(6,7) TITLE
    6 FORMAT(72A1)
    7 FORMAT(1X,72A1,/)
C.....READ IN THE INTITIAL CONCRETE TEMPERATURE AND THE CEQUIVALENT     
C.....CONCRETE SAND ROUGHNESS.  
      READ(5,*) TCONI,RSAND    
      WRITE(6,*) TCONI,RSAND 
C.....READ ICTC TO SPECIFY THE CONCRETE COMPOSITION.  ICTC=1 FOR
C.....DEFAULT LIMESTONE-COMMON SAND, ICTC=2 FOR DEFAULT SILICIOUS,      
C.....OR ICTC=3 FOR DEFAULT LIMESTONE-LIMESTONE.  IF ICTC=4, THE
C.....COMPOSITION IS USER-SPECIFIED.  IF ICTC=5, THE BASEMAT IS ALL 
C.....ALL STEEL.      
      READ(5,*) ICTC 
      WRITE(6,*) ICTC
      IF(ICTC.LE.3) GO TO 12    
      IF(ICTC.EQ.5) GO TO 122
C.....READ THE USER SPECIFIED COMPOSITIONS (WEIGHT PERCENT) IN THE      
C.....ORDER CO2, H2O, K2O, NA2O, TIO2, SIO2, CAO, MGO, AL2O3, FEO,      
C.....FE2O3, FE3O4, FE, AND CR. 
      READ(5,*) XWTC(1),XWTC(2),XWTC(3),XWTC(4),XWTC(5),XWTC(6),XWTC(7)
      WRITE(6,*) XWTC(1),XWTC(2),XWTC(3),XWTC(4),XWTC(5),XWTC(6),XWTC(7)
      READ(5,*) XWTC(8),XWTC(9),XWTC(10),XWTC(11),XWTC(12),XWTC(13),XWT
     1C(14) 
      WRITE(6,*) XWTC(8),XWTC(9),XWTC(10),XWTC(11),XWTC(12),XWTC(13),XWT
     1C(14)
C.....FOR THE USER-SPECIFIED CONCRETE COMPOSITION, READ THE CONCRETE    
C.....SOLIDUS, LIQUIDUS, AND DECOMPOSITION TEMPERATURES.
      READ(5,*) TCS,TCL,TDC
      WRITE(6,*) TCS,TCL,TDC
   12 CONTINUE
C.....FOR EITHER THE SPECIFIED OR DEFAULT CONCRETE TYPES, READ THE      
C.....TEMPERATURE RANGES FOR DRYOUT OF FREE WATER,  BOUND WATER
C.....(CA(OH)2), CALCIUM CARBONATE (CACO3), AND DOLOMITE(MGCA(CO3)2).   
      READ(5,*) TFWS,TFWL      
      READ(5,*) TBWS,TBWL      
      READ(5,*) TMCAS,TMCAL    
      READ(5,*) TCAS,TCAL      
      WRITE(6,*) TFWS,TFWL     
      WRITE(6,*) TBWS,TBWL     
      WRITE(6,*) TMCAS,TMCAL   
      WRITE(6,*) TCAS,TCAL     
  122 CONTINUE
C.....READ IN POUR DATA. FIRST READ IN THE ASSUMED OXIDE AND METAL      
C.....PHASE SOLIDUS AND LIQUIDUS TEMPERATURES.      
      READ(5,*) TFOS,TFOL      
      READ(5,*) TFMS,TFML      
      WRITE(6,*) TFOS,TFOL     
      WRITE(6,*) TFMS,TFML     
C.....READ IN TIME-DEPENDENT POUR RATE, COMPOSITION, TEMPERATURE, AND   
C.....DECAY HEAT DATA.
      READ(5,*) NPOURS
      WRITE(6,*) NPOURS
      DO 19 K=1,NPOURS
      READ(5,*) TST(K),TSTOP(K),AINTP(K),BINTP(K),ADEC(K),BDEC(K)      
      WRITE(6,*) TST(K),TSTOP(K),AINTP(K),BINTP(K),ADEC(K),BDEC(K)     
C.....READ THE TOTAL NUMBER OF CONSTITUENTS DRAINING AT THIS TIMESTEP   
C.....AND THE POUR RATE COEFFICIENTS FOR EACH CONSTITUENT.    
      READ(5,*) NISTP
      WRITE(6,*) NISTP
      IF(NISTP.EQ.0) GO TO 19   
      DO 20 ITP=1,NISTP
      READ(5,*) ITPA,APOUR(ITPA,K),BPOUR(ITPA,K)   
      WRITE(6,*) ITPA,APOUR(ITPA,K),BPOUR(ITPA,K)  
   20 CONTINUE
   19 CONTINUE
C.....READ MELT JET STREAM BREAKUP ANALYSIS MODEL INPUT DATA: NJET=0
C.....=> NO ANALYSIS, NJET=1 => THERMAL HOMOGENIZATION ANALYSIS, AND
C.....NET=2 => JET BREAKUP ANALYSIS.  NJETD=# OF POINTS IN RPV HOLE
C.....DIAMETER VS. TIME INTERP TABLE, AND ERPV=INITIAL ELEVATION OF 
C.....HOLE DRAIN SITE OVER BASEMAT.  ALSO, FOR THE CASE IN WHICH 
C.....NJET=2, NJETND=NUMBER OF NODES INTO WHICH THE PARTICLES FROM 
C.....JET FRAGMENTATION ARE ASSUMED TO SETTLE AND ACCUMULATE ON TOP
C.....OF THE CORE DEBRIS OR BASEMAT, WHICHEVER IS ENCOUNTERED.
C.....IF NJET < 2, SETTING OF NJETND DOES NOT MATTER.
      READ(5,*) NJET,NJETD,NJETND,ERPV,POROSBED
      WRITE(6,*) NJET,NJETD,NJETND,ERPV,POROSBED
      IF(NJET.EQ.0) GO TO 711
      DO 712 K=1,NJETD
      READ(5,*) TJETT(K),DJETT(K)
      WRITE(6,*) TJETT(K),DJETT(K)
  712 CONTINUE
      IF(NJET.EQ.1) GO TO 711
      DO 813 K=1,NJETND
      READ(5,*) NP
      WRITE(6,*) NP
      IFLGJ(NP)=1
  813 CONTINUE
  711 CONTINUE
C.....READ IN MELT CONSTITUENT OVERWRITE DATA IF A SIMULANT MATERIAL 
C.....IS TO BE USED IN THE SREADING ANALYSIS
      READ(5,*) NOVHT,COVS,COVL,DHSOV,XMOLOV,ROVS,ROVL
      WRITE(6,*) NOVHT,COVS,COVL,DHSOV,XMOLOV,ROVS,ROVL
      READ(5,*) NOVTK,TKOVS,TKOVL
      WRITE(6,*) NOVTK,TKOVS,TKOVL
      READ(5,*) NOVUM,VISOV
      WRITE(6,*) NOVUM,VISOV
      READ(5,*) NOVEM,EMOV
      WRITE(6,*) NOVEM,EMOV
      READ(5,*) NOVSIG,SIGOV
      WRITE(6,*) NOVSIG,SIGOV
C.....INITIALIZE SPREADING NODALIZATION.  FOR AUTOMATED MK I  
C.....NODALIZATION, SET NGEOM=1.  FOR USER-SPECIFIED INPUT, SET NGEOM=2.
      READ(5,*) NGEOM
      WRITE(6,*) NGEOM
      GO TO (31,32),NGEOM
   31 CONTINUE
C.....AUTOMATED NODALIZATION INITIALIZED BELOW.     
      READ(5,*) NSMP,NPED,NDOR,NDOOR,NSHL,NANULS   
      WRITE(6,*) NSMP,NPED,NDOR,NDOOR,NSHL,NANULS  
C.....READ REACTOR DIMENSIONS.  
      NUMNOD=NSMP+NPED+NDOR+NSHL+NANULS   
      READ(5,*) RSUMP,ELSMP,RPED,TPED,WDOOR,RSHELL 
      WRITE(6,*) RSUMP,ELSMP,RPED,TPED,WDOOR,RSHELL
C.....READ MATERIAL DISPERSION RADIUS BENEATH RPV.  
      READ(5,*) RCOMP
      WRITE(6,*) RCOMP
C.....READ SUMP COVER PLATE LOGIC CONTROL VARIABLE AND THE THICKNESS    
C.....OF THE COVER PLATE.
      READ(5,*) NSMPCV,THCKCV  
      WRITE(6,*) NSMPCV,THCKCV 
C.....USER-SPECIFIED SUMP BOILOVER OPTION PARAMETERS SET BELOW.  SET    
C.....NBOIL=0 TO BYPASS THIS OPTION, OR NBOIL=1 TO PERFORM.  IF NBOIL=1 
C.....SET TEBOIL=TEMPERATURE OF MELT IN SUMP WHEN BOILOVER OCCURS.  SET 
C.....TMBOIL=TIME AT WHICH BOILOVER OF THE SUMP VOLUME OCCURS.  SET     
C.....XBOIL VECTOR EQUAL TO THE MASS FRACTION OF CONSTITUENTS IN SUMP   
C.....AT THE TIME OF BOILOVER.  THE MASS OF CORIUM OVER THE COVER
C.....PLATES IS ADDED TO THE INTIAL INVENTORY IN THE SUMP.    
      READ(5,*) NBOIL,TMBOIL,TEBOIL,VFINT
      WRITE(6,*) NBOIL,TMBOIL,TEBOIL,VFINT
      IF(NBOIL.EQ.0) GO TO 21   
      DO 22 KVP=1,15  
      READ(5,*) XBLT(KVP)      
   22 CONTINUE
   21 CONTINUE
C.....READ NSWALL, WHICH CONTROLS THE SIDEWALL HEAT TRANSFER IN THE     
C.....THE DRYWELL ANNULUS.      
      READ(5,*) NSWALL
      WRITE(6,*) NSWALL
C.....READ SPREADING ANGLE OUTSIDE DOORWAY.
      READ(5,*) ANGFAN
      WRITE(6,*) ANGFAN
      GO TO 33
   32 CONTINUE
C.....MANUAL INPUT DONE BELOW.  
      NDOOR=1
      NSMPCV=0
      NSWALL=0
      ARC(1)=0.D0     
C.....READ NUMBER OF SPREADING NODES.     
      READ(5,*) NUMNOD
      WRITE(6,*) NUMNOD
      HASUM=0.D0      
      DO 37 K=1,NUMNOD
C.....READ THE NODE SIZE, NODE RADIAL LOCATION, ARC LENGTH, AND NODE    
C.....AREA. 
      KP1=K+1
      READ(5,*) IFLGA(K),DXNODE(K),RAD(K),ARC(KP1),AREA(K),ELEVAT(K),
     1HCP(K)   
      WRITE(6,*) IFLGA(K),DXNODE(K),RAD(K),ARC(KP1),AREA(K),ELEVAT(K),  
     1HCP(K)
      ELO(K)=ELEVAT(K)
      IF(IFLGA(K).EQ.1) HASUM=HASUM+AREA(K)
   37 CONTINUE
      DO 38 K=1,NUMNOD
      XFACMS(K)=0.D0  
      IF(IFLGA(K).EQ.1) XFACMS(K)=1.D0/HASUM
   38 CONTINUE
   33 CONTINUE
C.....READ IN DATA FOR VISCOSITY CALCULATION.  READ NVTPE=1 TO USE 
C.....ISHII-ZUBER VISCOSITY CORRRELATION, OR NVTPE=2 TO USE RAMCACCIOTTI 
C.....CORRELATION.  SET NSOLTP=1 TO USE LINEAR SOLID FRACTION VARIATION
C.....FOR OXIDE PHASE, OR SET NSOLTP=2 TO USE USER-SUPPLIED FUNCTION. SET
C.....NOSLF=NO. OF PTS. IN SOLID PHASE FUNCTION IF NSOLTP=2.  SET ALPMAX   
C.....MAXIMUM SOLID FRACTION, AND SET CRAMCON=CONSTANT IN RAMACCIIOTTI 
C.....CORRELATION IF THAT ONE IS USED. 
C.....RELOCATE.
      READ(5,*) NVTPE,NSOLTP,NSOLF,ALPMAX,CRAMCON
      WRITE(6,*) NVTPE,NSOLTP,NSOLF,ALPMAX,CRAMCON
      IF(NSOLTP.EQ.1) GO TO 8820
      DO 8810 K=1,NSOLF
      READ(5,*) TNORM(K),FRCSOL(K)
      WRITE(6,*) TNORM(K),FRCSOL(K)
 8810 CONTINUE
 8820 CONTINUE
C.....SETUP SUBSTRATE NODALIZATION. FIRST READ NUMBER OF VERTICAL NODES.
      READ(5,*) NMVER
      WRITE(6,*) NMVER
C.....READ INTERNODAL SPACING AND CELL INTERFACE SPACING.     
      DO 41 K=1,NMVER 
      READ(5,*) XBCN(K),XDCN(K)
      WRITE(6,*) XBCN(K),XDCN(K)
   41 CONTINUE
C.....ASSIGN BASEMAT NODAL SPACING BEFORE CONTINUING.
      IAS=ICTC
      IF(ICTC.GT.3) IAS=1
      DO 84 K=1,NUMNOD
      DO 84 L=1,NMVER 
      XBTW(L,K)=XBCN(K)
      XDIST(L,K)=XDCN(K)
      NTYPMT(L,K)=IAS 
      NL(L,K)=1
      IF(ICTC.EQ.5) NL(L,K)=2
   84 CONTINUE
C.....FLAG NODES WHICH ARE STEEL, AND SET NSIMST=0 IF STEEL NODES ARE
C.....STEEL, OR SET NSIMST=1 IF STEEL NODES ARE SIMULANT MATERIAL.
      READ(5,*) NSTEEL,NSIMST
      WRITE(6,*) NSTEEL,NSIMST
      IF(NSTEEL.EQ.0.AND.ICTC.LE.4) GO TO 999 
      IF(ICTC.EQ.5) GO TO 999
      DO 66 K=1,NSTEEL
      READ(5,*) IX,IY
      WRITE(6,*) IX,IY
      NL(IY,IX)=2     
   66 CONTINUE
  999 CONTINUE
C.....READ IN STEEL SIMULANT PROPERTIES IF STEEL IS PRESENT AND THE 
C.....USER WANTS TO OVERWRITE THE PROPERTY DATA AND USE SOME OTHER 
C.....MATERIAL.
      IF(NSTEEL.EQ.0.AND.ICTC.LE.4) GO TO 797
      IF(NSIMST.EQ.0) GO TO 797
C.....SIMULANT STEEL IF THIS POINT IS REACHED; READ DATA FOR 
C.....OVERWRITE.
      READ(5,*) TSIS,TSIL,CSIS,CSIL,DHSIL
      WRITE(6,*) TSIS,TSIL,CSIS,CSIL,DHSIL
      READ(5,*) RSIS,RSIL,TKSIS,TKSIL,ESIL
      WRITE(6,*) RSIS,RSIL,TKSIS,TKSIL,ESIL
      READ(5,*) XCSI(1),XCSI(2),XCSI(3),XCSI(4) 
      WRITE(6,*) XCSI(1),XCSI(2),XCSI(3),XCSI(4) 
      READ(5,*) XCSI(5),XCSI(6),XCSI(7),XCSI(8) 
      WRITE(6,*) XCSI(5),XCSI(6),XCSI(7),XCSI(8) 
      READ(5,*) XCSI(9),XCSI(10),XCSI(11),XCSI(12) 
      WRITE(6,*) XCSI(9),XCSI(10),XCSI(11),XCSI(12) 
      READ(5,*) XCSI(13),XCSI(14),XCSI(15),XCSI(16) 
      WRITE(6,*) XCSI(13),XCSI(14),XCSI(15),XCSI(16) 
  797 CONTINUE
C.....SET NSKIPE=0 TO BYPASS SHELL LINER HEATUP CALCULATION.  
      READ(5,*) NSKIPE
      WRITE(6,*) NSKIPE
      IF(NSKIPE.EQ.0) GO TO 51  
C.....READ SHELL NODALIZATION & CONTROL DATA.
      READ(5,*) NBMADJ,NUMSHH,NUMSHV,NBFZOE,NCRTEM,NLOGSH    
      WRITE(6,*) NBMADJ,NUMSHH,NUMSHV,NBFZOE,NCRTEM,NLOGSH
      NBFZME=NBFZOE   
C.....READ SHELL INITIAL TEMPERATURE, HEIGHT TO BE NODALIZED, SHELL     
C.....THICKNESS, SHELL ANGLE WITH RESPECT TO VERTICAL, AND WIDTH OF    
C....."SLOT JET" IMPINGING ON SHELL.
      READ(5,*) TSHELI,HNODT,THSHL,ANGSHL,BWIDTH  
      WRITE(6,*) TSHELI,HNODT,THSHL,ANGSHL,BWIDTH 
   51 CONTINUE
C.....READ OVERLYING MEDIUM HEAT TRANSFER PROPERTIES.  SET NWAT=0 IF    
C.....WATER IS ABSENT FROM THE CAVITY, SET NWAT=1 IF WATER IS PRESENT
C.....OVER BASEMAT AT A CONSTANT TEMP AND LEVEL LIMITED BY DOWNCOMER
C.....HEIGHT, OR SET NWAT=2 IF DETAILED WATER INVENTORY MODEL IS USED.
C.....IF NWAT=2 THEN TWATI AND ELWATI ARE INITIAL WATER TEMPERATURE
C.....AND DEPTH, RESPECTIVELY. 
      READ(5,*) TBOUND,EMISCN,PDRYWL     
      WRITE(6,*) TBOUND,EMISCN,PDRYWL   
      READ(5,*) NWAT,HDOWNC,TWATI,ELWATI
      WRITE(6,*) NWAT,HDOWNC,TWATI,ELWATI
C.....IF NWAT<2, THEN BYPASS THE WATER INJECTION AND SPILLOVER INPUT
C.....DATA.  IF NWAT=2, READ IN INJECTION AND SPILLOVER LOCATION/SIZE
C.....DATA.
      IF(NWAT.LT.2) GO TO 199
      READ(5,*) NINJ
      WRITE(6,*) NINJ
C.....READ WATER INJECTION DATA: NINJ=# OF INJECTION POINTS, NINJP ARE 
C.....NODES WHERE INJECTION IS OCCURING, NPTINJ IS # OF POINTS IN TEMP 
C.....& FLOWRATE INTERPOLATION TABLE AT EACH INJECTION POINT, AND  
C.....TIMINJ, XDTINJ, AND TDTINJ ARE TIME, FLOWRATE, AND TEMP POINTS. 
      IF(NINJ.EQ.0) GO TO 189
      DO 219 K=1,NINJ
      READ(5,*) NINJP(K),NPTINJ(K)
      WRITE(6,*) NINJP(K),NPTINJ(K)
      DO 220 L=1,NPTINJ(K)
      READ(5,*) TIMINJ(K,L),XDTINJ(K,L),TDTINJ(K,L)
      WRITE(6,*) TIMINJ(K,L),XDTINJ(K,L),TDTINJ(K,L)
  220 CONTINUE
  219 CONTINUE
  189 CONTINUE
      READ(5,*) NDOWNC
      WRITE(6,*) NDOWNC
      IF(NDOWNC.EQ.0) GO TO 199
C.....READ SPILLOVER DATA: NDOWNC=# OF SPILLOVER POINTS, NDC ARE 
C.....NODES WHERE SPILLOVER IS OCCURING, ELDCO IS INLET ELEVATION AT THE 
C.....LOCATION, AND NPTDC IS # OF INTERPOLATION POINTS IN THE SPILLOVER 
C.....CORD LENGTH TABLE AS A FUNCTION OF ELEVATION (ELDCX,CRDCX).
      DO 319 K=1,NDOWNC
      READ(5,*) NDC(K),ELDCO(K),NPTDC(K)
      WRITE(6,*) NDC(K),ELDCO(K),NPTDC(K)
      DO 320 L=1,NPTDC(K)
      READ(5,*) ELDCX(K,L),CRDCX(K,L)
      WRITE(6,*) ELDCX(K,L),CRDCX(K,L)
  320 CONTINUE
  319 CONTINUE
  199 CONTINUE      
C.....READ IN CRUST, ABLATION, & HEAT TRANSFER COEFFICIENT OPTIONS.     
      READ(5,*) NBFRZO,NTHINC,NABLFM,NDRNFM,NBCBOT,NCRTOP,NADAB    
      WRITE(6,*) NBFRZO,NTHINC,NABLFM,NDRNFM,NBCBOT,NCRTOP,NADAB   
      IF(NADAB.EQ.1) NBFRZO=1
      NBFRZM=NBFRZO   
C.....READ MELT/BASEMAT INTERFACIAL HEAT TRANSFER RESISTANCE CONTROL
C.....PARAMETER NINTF=0 IF NO INTERFACIAL HEAT TRANSFER RESISTANCE,
C.....OR NINTF=1 IF RESISTANCE IS HINTF.
      READ(5,*) NINTF,HINTF
      WRITE(6,*) NINTF,HINTF
C.....READ ABLATION MAPPING CRITERIA, AND DESIRED NODE DEPTH FOR
C.....SOLIDIFIED DEBRIS.
      READ(5,*) XFCABL,XNDMIN  
      WRITE(6,*) XFCABL,XNDMIN 
C.....READ INTEGRATION CONTROL INFORMATION FOR MELT FLUID MECHANICS
C.....CALCULATION.    
      READ(5,*) NVELP,NITMAX,DAVMX,DVMX,NINVIS
      WRITE(6,*) NVELP,NITMAX,DAVMX,DVMX,NINVIS
C.....READ MELT SPECIFIC ENTHALPY INTEGRATION CONTROL DATA    
      READ(5,*) NENMAX,DEAVMX,DEMX
      WRITE(6,*) NENMAX,DEAVMX,DEMX      
C.....READ INTEGRATION CONTROL INFORMATION FOR WATER FLUID MECHANICS
C.....CALCULATION (SETTINGS ARE ARBITRARY IF NWAT<2).
      READ(5,*) NVELPW,NITMAXW,DAVMXW,DVMXW
      WRITE(6,*) NVELPW,NITMAXW,DAVMXW,DVMXW
C.....READ WATER SPECIFIC ENTHALPY INTEGRATION CONTROL DATA    
      READ(5,*) NENMXW,DEAVMXW,DEMXW
      WRITE(6,*) NENMXW,DEAVMXW,DEMXW      
C.....READ PRINTOUT AND PLOTTING CONTROL DATA: NPRINT=1 FOR NORMAL
C.....OUTPUT, OR NPRINT=2 FOR VERBOSE; NPFREQ=HOW MANY TIMESTEPS
C.....BETWEEN ASCII OUTPUT, AND NPEND=HOW MANY SHELL NODES 
C.....HORIZONTALLY TO PRINT, IF APPLICABLE.  ALSO SET NBEDCQ=0 IF
C.....PARTICULATE FORMED FROM JET EROSION (NJET=2) IS TO BE 
C.....RETAINED AS A SEPARATE DEBRIS BED IN CORQUENCH INITIALIZATION 
C.....FILE, OR NBEDCQ=1 IF IT IS TO BE INTEGRATED WITH UNDERLYING
C.....MELT AND THERMALLY EQUILIBRATED WITH THAT MATERIAL.  THE
C.....SETTING OF THIS PARAMETR IS NOT IMPORTANT IF THE OUTPUT
C.....DATA IN CQINPUT DOES NOT MATTER.
      READ(5,*) NPRINT,NPFREQ,NPEND,NBEDCQ
      WRITE(6,*) NPRINT,NPFREQ,NPEND,NBEDCQ
C.....READ HOW MANY TIMES THE SPATIALLY-DEPENDENT DATA ARE TO BE
C.....WRITTEN OUT AND AT WHAT TIMES THE DATA ARE TO BE WRITTEN.
      READ(5,*) NTIMSPC
      WRITE(6,*) NTIMSPC
      IF(NTIMSPC.EQ.0) GO TO 733
      DO 734 KP=1,NTIMSPC
      READ(5,*) TIMSPC(KP)
      WRITE(6,*) TIMSPC(KP)
  734 CONTINUE
  733 CONTINUE
C.....READ AT WHAT FREQUENCY THE DATA ARE TO BE PRINTED OUT AT
C.....SELECTED NODES, THE TOTAL NUMBER OF NODES TO BE PRINTED,
C.....AND THE SPECFIC NODES AT WHICH THE DATA ARE TO BE PRINTED.
      READ(5,*) NPLFREQ,NPLTOT
      WRITE(6,*) NPLFREQ,NPLTOT 
      IF(NPLTOT.EQ.0) GO TO 333
      DO 334 KP=1,NPLTOT
      READ(5,*) NPLLOC(KP)
      WRITE(6,*) NPLLOC(KP)
  334 CONTINUE
  333 CONTINUE
C.....READ THE NUMBER OF BASEMAT NODES WHERE THERMAL RESPONSE IS TO
C.....BE PLOTTED AND THE NODE LOCATIONS (HORIZONTAL POSITION, 
C.....VERTICAL POSITION) 
      READ(5,*) NBPL
      WRITE(6,*) NBPL
      IF(NBPL.EQ.0) GO TO 453
      DO 454 KD=1,NBPL
      READ(5,*) IXP(KD),IYP(KD)
  454 CONTINUE
  453 CONTINUE
C.....READ THE INITIAL TIME, INTEGRATION TIMESTEP, AND MAXIMUM TIME     
C.....TO WHICH THE INTEGRATION IS TO BE PERFORMED.  
      READ(5,*) TIMEO,DTIME,TMAX
      WRITE(6,*) TIMEO,DTIME,TMAX
      WRITE(6,951)
  951 FORMAT(//,1X,'END ECHO-WRITE OF INPUT; START CALCULATION',//)
      RETURN
      END