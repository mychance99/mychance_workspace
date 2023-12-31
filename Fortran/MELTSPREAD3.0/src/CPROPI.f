C...................................................................... 
C.....THIS ROUTINE DEFINES VARIABLES REQUIRED TO DETERMINE THE
C.....CONCRETE TEMPERATURE GIVEN THE CONCRETE SPECIFIC ENTHALPY.
C...................................................................... 
      SUBROUTINE CPROPI
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
      COMMON/DATSH/TSTLS,TSTLL,ESTLS,ESTLL,DFUSST,CSTLL
      COMMON/DATSIMS/TSIS,TSIL,CSIS,CSIL,DHSIL,RSIS,RSIL,TKSIS,TKSIL,
     1ESIL,XCSI(16)
      COMMON/DATOVER/COVS,COVL,DHSOV,XMOLOV,ROVS,ROVL,TKOVS,TKOVL,VISOV,
     1EMOV,SIGOV
      DATA ZERO / 0.0 /
      DATA ONE / 1.0 /
      DATA TWO / 2.0 /
      DATA HALF / 0.5 /
      DATA CAL / 4.184 /
      DATA GRAM / 0.001 /
      DATA TREF/298.15D0/
C.....SET DATA AND INTEGER FLAGS.
      ILCS=1
      IBAS=2
      ILL=3 
      GRAV=9.82D0     
      STEF=5.67D-8    
      PI=3.141592654D0
      XMH2O=18.0152D0*GRAM      
      XMCO2=44.0098D0*GRAM      
      XMCACO=100.0892D0*GRAM    
      XMMGCA=184.4034D0*GRAM    
      XMCAOH=74.0946D0*GRAM     
      QFEH2O=.078D6   
      QCRH2O=3.57D6   
      QZRH2O=6.74D6   
      QFECO2=0.D0     
      QCRCO2=2.75D6   
      QZRCO2=5.84D6   
      XZRMRE=2.D0     
      XFEMRE=1.D0     
      XCRMRE=1.5D0    
      XZRORE=2.D0     
      XFEORE=1.D0     
      XCRORE=3.D0     
C.....ASSIGN THE MATERIAL INDICES.
      NMATC=18
      ICAOH2=1
      ICACO3=2
      IMCCO3=3
      IFH2O=4
      IVH2O=5
      ICK2O=6
      IVK2O=7
      INA2O=8
      ITIO2=9
      ISIO2=10
      ICAO=11
      IMGO=12
      IAL2O3=13
      IFEO=14
      IFE2O3=15
      IFE3O4=16
      IFE=17
      ICR=18
      INI=19
      IZR=20
      IU=21 
      IB4C=22
      IB=23 
      ICR2O3=24
      INIO=25
      IUO2=26
      IZRO2=27
      IB2O3=28
C.....DEFINE THE MOLECULAR WEIGHT, ROOM TEMPERATURE DENSITY,  
C.....SOLID PHASE SPECIFIC ENTHALPY EQUATION COEFFICIENTS, AND
C.....LIQUID PHASE SPECIFIC ENTHALPY EQUATION COEFFICIENTS FOR
C.....EACH MATERIAL.  
      I=ISIO2
      IMOX(I)=1
      FMMOL(I)=60.0848*GRAM     
      ROM(I)=2.20*1000.0
      ROMLIQ(I)=2.129537062*1000.0
      AEQM(I,1)=1.269798875D-03*CAL
      BEQM(I,1)=13.84284222*CAL 
      CEQM(I,1)=-4.2379302D+03*CAL
      AEQM(I,2)=ZERO  
      BEQM(I,2)=20.5*CAL
      CEQM(I,2)=-10096.00*CAL   
      I=ICAO
      IMOX(I)=1
      FMMOL(I)=56.0794*GRAM     
      ROM(I)=3.25*1000.0
      ROMLIQ(I)=2.878868457*1000.0
      AEQM(I,1)=6.7956134D-04*CAL
      BEQM(I,1)=11.25975464*CAL 
      CEQM(I,1)=-3.415754648D+03*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=15*CAL
      CEQM(I,2)=10352.0*CAL     
      I=IMGO
      IMOX(I)=1
      FMMOL(I)=40.3114*GRAM     
      ROM(I)=3.58*1000.0
      ROMLIQ(I)=3.143761653*1000.0
      AEQM(I,1)=6.054886725D-04*CAL
      BEQM(I,1)=10.71904531*CAL 
      CEQM(I,1)=-3.248045318D+03*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=14.5*CAL
      CEQM(I,2)=9099.0*CAL      
      I=ICAOH2
      IMOX(I)=1
      FMMOL(I)=74.09474*GRAM    
      ROM(I)=2.24*1000.0
      ROMLIQ(I)=ROM(I)
      AEQM(I,1)=2.95601497D-03*CAL
      BEQM(I,1)=21.65597006*CAL 
C.....CSOL COEFFICIENT INCORRECT FOR ICAOH2; CORRECTED 11/25/92.
      CEQM(I,1)=-32875.98*CAL
      AEQM(I,2)=ZERO  
      BEQM(I,2)=ZERO  
      CEQM(I,2)=ZERO  
      I=ICACO3
      IMOX(I)=1
      FMMOL(I)=100.08935*GRAM   
      ROM(I)=2.710*1000.0
      ROMLIQ(I)=ROM(I)
      AEQM(I,1)=4.012976829D-03*CAL
      BEQM(I,1)=21.20585561*CAL 
C.....CSOL COEFFICIENT INCORRECT FOR ICACO3; CORRECTED 11/25/92.
      CEQM(I,1)=-4.91757D+04*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=ZERO  
      CEQM(I,2)=ZERO  
      I=IMCCO3
      IMOX(I)=1
      FMMOL(I)=184.4107*GRAM    
      ROM(I)=2.872*1000.0
      ROMLIQ(I)=ROM(I)
      AEQM(I,1)=1.082974027D-02*CAL
      BEQM(I,1)=34.06862336*CAL 
      CEQM(I,1)=-1.77764174D+05*CAL
      AEQM(I,2)=ZERO  
      BEQM(I,2)=ZERO  
      CEQM(I,2)=ZERO  
      I=IFH2O
      IMOX(I)=1
      FMMOL(I)=18.01534*GRAM    
      ROM(I)=0.9970636476*1000.0
      ROMLIQ(I)=ROM(I)
      AEQM(I,1)=1.030120773D-03*1000.0*FMMOL(I)     
      BEQM(I,1)=3.485859882*1000.0*FMMOL(I)
C.....CSOL COEFFICIENT INCORRECT FOR IFH2O; CORRECTED 11/25/92.
      CEQM(I,1)=-1.113026D+03*1000.0*FMMOL(I)    
      AEQM(I,2)=ZERO  
      BEQM(I,2)=ZERO  
      CEQM(I,2)=ZERO  
      I=IVH2O
      IMOX(I)=1
      FMMOL(I)=ZERO   
      ROM(I)=ROM(IFH2O)
      ROMLIQ(I)=ROM(I)
      AHVAP=-4.737086585D-03*1000.0*FMMOL(IFH2O)    
      BHVAP=7.57401951D-01*1000.0*FMMOL(IFH2O)      
      CHVAP=2.637717168D+03*1000.0*FMMOL(IFH2O)     
      AEQM(I,1)=ZERO  
      BEQM(I,1)=ZERO  
      CEQM(I,1)=AEQM(IFH2O,1)*TFWL*TFWL+BEQM(IFH2O,1)*TFWL+CEQM(IFH2O,1)
     1+AHVAP*TFWL*TFWL+BHVAP*TFWL+CHVAP   
      AEQM(I,2)=ZERO  
      BEQM(I,2)=ZERO  
      CEQM(I,2)=CEQM(I,1)
      I=IAL2O3
      IMOX(I)=1
      FMMOL(I)=101.9612*GRAM    
      ROM(I)=3.965*1000.0
      ROMLIQ(I)=3.741451713*1000.0
      AEQM(I,1)=6.587394436D-03 
      BEQM(I,1)=108.2762663     
      CEQM(I,1)=-3.286814483D+04
      AEQM(I,2)=ZERO  
      BEQM(I,2)=192.464
      CEQM(I,2)=-82016.728      
      I=IFEO
      IMOX(I)=1
      FMMOL(I)=71.8464*GRAM     
      ROM(I)=5.7*1000.0
      ROMLIQ(I)=5.32912083*1000.0
      AEQM(I,1)=1.1278054D-03*CAL
      BEQM(I,1)=11.78824218*CAL 
      CEQM(I,1)=-3.6130498D+03*CAL
      AEQM(I,2)=ZERO  
      BEQM(I,2)=16.3*CAL
      CEQM(I,2)=-2237.0*CAL     
      I=IFE2O3
      IMOX(I)=1
      FMMOL(I)=159.6922*GRAM    
      ROM(I)=5.24*1000.0
      ROMLIQ(I)=4.949536681*1000.0
      AEQM(I,1)=3.452197061D-04*CAL
      BEQM(I,1)=33.56608762*CAL 
      CEQM(I,1)=-29133.351*CAL  
      AEQM(I,2)=ZERO  
      BEQM(I,2)=34.0*CAL
      CEQM(I,2)=9356.0*CAL      
      I=INA2O
      IMOX(I)=1
      FMMOL(I)=61.979*GRAM      
      ROM(I)=2.27*1000.0
      ROMLIQ(I)=ROM(I)
      AEQM(I,1)=-5.86089669D-04*CAL
      BEQM(I,1)=25.45791197*CAL 
      CEQM(I,1)=-7.53441066D+03*CAL
      AEQM(I,2)=ZERO  
      BEQM(I,2)=25.0*CAL
      CEQM(I,2)=3352.0*CAL      
      I=ICK2O
      IMOX(I)=1
      FMMOL(I)=94.2034*GRAM     
      ROM(I)=2.32*1000.0
      ROMLIQ(I)=ROM(I)
      AEQM(I,1)=5.003493755D-03*CAL
      BEQM(I,1)=17.24193641*CAL 
      CEQM(I,1)=-9.238242731D+04*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=ZERO  
      CEQM(I,2)=ZERO  
      I=IVK2O
      IMOX(I)=1
      FMMOL(I)=ZERO   
      ROM(I)=ROM(ICK2O)
      ROMLIQ(I)=ROM(I)
      AEQM(I,1)=ZERO  
      BEQM(I,1)=ZERO  
      CEQM(I,1)=ZERO  
      AEQM(I,2)=ZERO  
      BEQM(I,2)=ZERO  
      CEQM(I,2)=CEQM(ICK2O,1)+TCAL*(BEQM(ICK2O,1)+TCAL*AEQM(ICK2O,1))   
     1+37024.8566*CAL 
      I=ITIO2
      IMOX(I)=1
      FMMOL(I)=79.8988*GRAM     
      ROM(I)=4.26*1000.0
      ROMLIQ(I)=4.020841193*1000.0
      AEQM(I,1)=2.717132431D-03 
      BEQM(I,1)=67.65301585     
      CEQM(I,1)=-2.041228188D+04
      AEQM(I,2)=ZERO  
      BEQM(I,2)=100.416
      CEQM(I,2)=-10926.08
      I=IFE 
      IMOX(I)=0
      FMMOL(I)=55.847*GRAM      
      ROM(I)=7.867*1000.0
      ROMLIQ(I)=7.01*1000.0     
      AEQM(I,1)=5.928784491D-04*CAL
      BEQM(I,1)=8.024965771*CAL 
      CEQM(I,1)=-2.445346539D+03*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=11.0*CAL
      CEQM(I,2)=-2587.0*CAL     
      I=ICR 
      IMOX(I)=0
      FMMOL(I)=51.996*GRAM      
      ROM(I)=7.19*1000.0
      ROMLIQ(I)=6.28*1000.0     
      AEQM(I,1)=1.839387325D-03*CAL
      BEQM(I,1)=4.264209995*CAL 
      CEQM(I,1)=-1.434883645D+03*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=9.4*CAL
      CEQM(I,2)=18.0*CAL
      I=INI 
      IMOX(I)=0
      FMMOL(I)=58.71*GRAM
      ROM(I)=8.908*1000.0
      ROMLIQ(I)=7.77*1000.0     
      AEQM(I,1)=4.865234897D-04*CAL
      BEQM(I,1)=6.970520914*CAL 
      CEQM(I,1)=-2.121509549D+03*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=10.3*CAL
      CEQM(I,2)=-2242.8*CAL     
      I=IFE3O4
      IMOX(I)=1
      FMMOL(I)=231.5386*GRAM    
      ROM(I)=5.18*1000.0
      ROMLIQ(I)=4.725017173*1000.0
      AEQM(I,1)=ZERO  
      BEQM(I,1)=48.60305344*CAL 
      CEQM(I,1)=-14483.70992*CAL
      AEQM(I,2)=ZERO  
      BEQM(I,2)=51.0*CAL
      CEQM(I,2)=14034.0*CAL     
      I=IZR 
      IMOX(I)=0
      FMMOL(I)=91.22*GRAM
      ROM(I)=6.50*1000.0
      ROMLIQ(I)=6.06*1000.0     
      AEQM(I,1)=ZERO  
      BEQM(I,1)=7.903221392*CAL 
      CEQM(I,1)=-2.356345458D+03*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=8.0*CAL
      CEQM(I,2)=1476.0*CAL      
      I=IU  
      IMOX(I)=0
      FMMOL(I)=238.03*GRAM      
      ROM(I)=19.0*1000.0
      ROMLIQ(I)=17.905*1000.0   
      AEQM(I,1)=ZERO  
      BEQM(I,1)=10.59764196*CAL 
      CEQM(I,1)=-3.15968695D+03*CAL
      AEQM(I,2)=ZERO  
      BEQM(I,2)=11.45*CAL
      CEQM(I,2)=-2.32125D+03*CAL
      I=IB4C
      IMOX(I)=0
      FMMOL(I)=55.25515*GRAM    
      ROM(I)=2.52*1000.0
      ROMLIQ(I)=2.395891505*1000.0
      AEQM(I,1)=1.327567331D-02 
      BEQM(I,1)=84.7386562      
      CEQM(I,1)=-2.644495038D+04
      AEQM(I,2)=ZERO  
      BEQM(I,2)=135.98
      CEQM(I,2)=37486.86
      I=IB  
      IMOX(I)=0
      FMMOL(I)=10.811*GRAM      
      ROM(I)=2.37*1000.0
      ROMLIQ(I)=2.08*1000.0     
      AEQM(I,1)=8.148631809D-04*CAL
      BEQM(I,1)=4.081629368*CAL 
      CEQM(I,1)=-1.289373773D+03*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=7.5*CAL
      CEQM(I,2)=7159.0*CAL      
      I=ICR2O3
      IMOX(I)=1
      FMMOL(I)=151.9902*GRAM    
      ROM(I)=5.21*1000.0
      ROMLIQ(I)=4.941650398*1000.0
      AEQM(I,1)=5.467302068D-03 
      BEQM(I,1)=115.0902254     
      CEQM(I,1)=-3.480015791D+04
      AEQM(I,2)=ZERO  
      BEQM(I,2)=156.9 
      CEQM(I,2)=23117.3
      I=INIO
      IMOX(I)=1
      FMMOL(I)=74.7094*GRAM     
      ROM(I)=6.67*1000.0
      ROMLIQ(I)=6.043098796*1000.0
      AEQM(I,1)=9.374674435D-04*CAL
      BEQM(I,1)=11.50727196*CAL 
      CEQM(I,1)=-3.512417903D+03*CAL      
      AEQM(I,2)=ZERO  
      BEQM(I,2)=14.30956023*CAL 
      CEQM(I,2)=7046.353151*CAL 
      I=IUO2
      IMOX(I)=1
      FMMOL(I)=270.0288*GRAM    
      ROM(I)=10.97*1000.0
      ROMLIQ(I)=8.739*1000.0    
      AEQM(I,1)=1.816681331D-02 
      BEQM(I,1)=40.90328306     
      CEQM(I,1)=-8.077459866D+03
      AEQM(I,2)=ZERO  
      BEQM(I,2)=130.95
      CEQM(I,2)=-30911.0
      I=IZRO2
      IMOX(I)=1
      FMMOL(I)=123.2188*GRAM    
C$$$$$ZrO2 SOLIDUS/LIQUIDUS DENSITIES SWITCHED ON 1-8-96 BY MTF.
C$$$$$ZrO2 DENSITY AT LIQUIDUS IS ACTUALLY GREATER THAN SOLIDUS.
      ROM(I)=5.8*1000.0      
      ROMLIQ(I)=5.9914*1000.0      
      AEQM(I,1)=ZERO  
      BEQM(I,1)=0.6937265932D+03*FMMOL(I) 
      CEQM(I,1)=-2.068345838D+05*FMMOL(I) 
      AEQM(I,2)=ZERO  
      BEQM(I,2)=0.815D+03*FMMOL(I)
      CEQM(I,2)=1.39D+05*FMMOL(I)
C.....INPUT B2O3 DATA, UNLESS A SIMULANT MATERIAL IS BEING USED,
C.....IN WHICH CASE OVERWRITE THIS DATA.
      I=IB2O3
      IMOX(I)=1
      IF(NOVHT.EQ.0) GO TO 6443
C.....PROPERTY DATA REASSIGNED HERE.
      FMMOL(I)=XMOLOV*GRAM     
      ROM(I)=ROVS
      ROMLIQ(I)=ROVL
      AEQM(I,1)=ZERO
      BEQM(I,1)=COVS*FMMOL(I)
      CEQM(I,1)=-TREF*COVS*FMMOL(I)
      AEQM(I,2)=ZERO  
      BEQM(I,2)=COVL*FMMOL(I)
      CEQM(I,2)=(COVS*(TFOS-TREF)+DHSOV-COVL*TFOL)*FMMOL(I)
      GO TO 6444
C.....ACTUAL B2O3 DATA DONE HERE.
 6443 CONTINUE 
      FMMOL(I)=69.6202*GRAM     
      ROM(I)=1.812*1000.0
      ROMLIQ(I)=ROM(I)
      AEQM(I,1)=4.288754823D-02 
      BEQM(I,1)=44.67660527     
      CEQM(I,1)=-1.71327508D+04 
      AEQM(I,2)=ZERO  
      BEQM(I,2)=129.704
      CEQM(I,2)=-32118.922      
 6444 CONTINUE
C.....ASSIGN PROPERTIES FOR CARBON STEEL OR SIMULANT, 
C.....AS APPROPRIATE. 
      IF(NSTEEL.EQ.0.AND.ICTC.LE.4) GO TO 2400 
      IF(NSIMST.EQ.0) GO TO 2400
      TSTLS=TSIS
      TSTLL=TSIL
      ESTLS=CSIS*(TSTLS-298.15) 
      ESTLL=ESTLS+DHSIL
      DFUSST=DHSIL
      CSTLL=CSIL
      GO TO 2401
 2400 CONTINUE
      TSTLS=1810.D0   
      TSTLL=1811.D0   
      ESTLS=1.0454D6  
      ESTLL=1.293D6   
      DFUSST=2.476D5  
      CSTLL=835.D0    
 2401 CONTINUE
C.....ASSIGN CONCRETE AND SLAG W/O'S, AND CONCRETE MOLE %'S, FOR THE    
C.....VARIOUS TYPES OF CONCRETE.
      DO 5807 KTP=1,NMATC
      XMOL(KTP)=0.D0  
 5807 CONTINUE
      IF(ICTC.EQ.4) GO TO 1135  
      ICTYPE=ICTC     
      IF(ICTC.EQ.5) ICTYPE=1
      GO TO(1131,1132,1133),ICTYPE
 1131 CONTINUE
      XFRGAS=.28     
      TCS=1393.D0     
      TCL=1568.D0     
      WPCC=.3425      
      WPM=0.1164      
      WPA=.0359
      WPCS=0.5050     
      XFH2OU=0.0429  
      XFCAOH=0.0840   
      XFMGCA=0.4579   
      XFCACO=.0093    
      XMOL(ISIO2)=.3821
      XMOL(IMGO)=.0000
      XMOL(ICAO)=.0857
      XMOL(IFEO)=.0000
      XMOL(IAL2O3)=.0281
      XMOL(INA2O)=.0141
      XMOL(ICK2O)=0.0051
      XMOL(ITIO2)=.0014
      XMOL(ICAOH2)=.0885
      XMOL(ICACO3)=.0073
      XMOL(IMCCO3)=.1938
      XMOL(IFH2O)=.1858
      XMOL(IVH2O)=0.
      GO TO 1134      
 1132 CONTINUE
      XFRGAS=.0791
      TCS=1403.D0     
      TCL=1523.D0     
      WPCC=0.1227     
      WPM=0.0057
      WPA=0.0276
      WPCS=.8440      
      XFH2OU=.0504  
      XFCAOH=0.0844   
      XFMGCA=.0215    
      XFCACO=0.0000  
      XMOL(ISIO2)=.6574
      XMOL(ICAO)=.0691
      XMOL(IAL2O3)=.0222
      XMOL(ICK2O)=.0084
      XMOL(IFE2O3)=.0035
      XMOL(ITIO2)=.0057
      XMOL(IMGO)=.0034
      XMOL(INA2O)=.0064
      XMOL(ICAOH2)=.0629
      XMOL(IMCCO3)=.0064
      XMOL(IFH2O)=.0154
      XMOL(IVH2O)=.00
      GO TO 1134      
 1133 CONTINUE
      XFRGAS=.406
      TCS=1495.D0     
      TCL=2577.D0     
      WPCC=0.7066     
      WPM=0.1114      
      WPA=0.0244
      WPCS=0.1577     
      XFH2OU=0.0503   
      XFCAOH=0.0843   
      XFMGCA=.3424    
      XFCACO=.4024   
      XMOL(INA2O)=.0000
      XMOL(ICK2O)=.00362
      XMOL(IMGO)=.0000
      XMOL(ICAO)=.0611
      XMOL(IAL2O3)=.0159
      XMOL(ISIO2)=.0993
      XMOL(ITIO2)=.00107
      XMOL(IFE2O3)=.0043
      XMOL(ICAOH2)=.0946
      XMOL(ICACO3)=.3341
      XMOL(IMCCO3)=.1542
      XMOL(IVH2O)=.2318
      XMOL(IFH2O)=.0
      GO TO 1134      
 1135 CONTINUE
      ICTYPE=1
C.....USER-SPECIFIED CONCRETE INTIALIZED BELOW. FIRST FIND SLAG
C.....DENSITY AND VOLUME FRACTIONS OF SLAG CONSTITUENTS USING THE
C.....BRUTE FORCE TECHNIQUE.    
      DENOM=XWTC(6)/FMMOL(ISIO2)+XWTC(7)/FMMOL(ICAO)+XWTC(8)/FMMOL(IMGO)
     1+XWTC(9)/FMMOL(IAL2O3)    
      FMSI=XWTC(6)/(FMMOL(ISIO2)*DENOM)   
      FMCA=XWTC(7)/(FMMOL(ICAO)*DENOM)    
      FMMG=XWTC(8)/(FMMOL(IMGO)*DENOM)    
      FMAL=XWTC(9)/(FMMOL(IAL2O3)*DENOM)  
      DENOM2=(FMSI*FMMOL(ISIO2))/ROMLIQ(ISIO2)      
     1+(FMCA*FMMOL(ICAO))/ROMLIQ(ICAO)    
     2+(FMMG*FMMOL(IMGO))/ROMLIQ(IMGO)    
     3+(FMAL*FMMOL(IAL2O3))/ROMLIQ(IAL2O3)
      WPCS=(FMSI*FMMOL(ISIO2))/(ROMLIQ(ISIO2)*DENOM2)
      WPA=(FMAL*FMMOL(IAL2O3))/(ROMLIQ(IAL2O3)*DENOM2)
      WPM=(FMMG*FMMOL(IMGO))/(ROMLIQ(IMGO)*DENOM2)  
      WPCC=(FMCA*FMMOL(ICAO))/(ROMLIQ(ICAO)*DENOM2) 
C.....SET THE TOTAL GAS FRACTION (WT%) IN THE INITIAL CONCRETE.
      XFRGAS=0.01D0*(XWTC(1)+XWTC(2))
C.....DECOMPOSE CO2 & H2O INTO FREE WATER, CAOH2, CACO3, AND DOLOMITE.  
      XHSAVE=XWTC(2)  
      HSUMT=ZERO      
      DO 61 K=1,14    
      IF(K.EQ.1) XML=XMCO2      
      IF(K.EQ.2) XML=FMMOL(K+2) 
      IF(K.EQ.3) XML=FMMOL(K+3) 
      IF(K.GE.4) XML=FMMOL(K+4) 
      HSUMT=HSUMT+XWTC(K)/XML   
   61 CONTINUE
      DO 62 K=1,14    
      IF(K.EQ.1) XML=XMCO2      
      IF(K.EQ.2) XML=FMMOL(K+2) 
      IF(K.EQ.3) XML=FMMOL(K+3) 
      IF(K.GE.4) XML=FMMOL(K+4) 
      XWTC(K)=XWTC(K)/(HSUMT*XML)
   62 CONTINUE
C.....PUT H2O INTO FREE WATER AND CAOH2; BOUND WATER (CAOH2) IS ASSUMED
C.....TO AMOUNT TO 2 WT% OF THE H2O INVENTORY.
      IF(XWTC(2).EQ.ZERO) GO TO 63
      FRACF=TWO/XHSAVE
      IF(FRACF.GT.ONE) FRACF=ONE
      FARG=XWTC(7)/XWTC(2)
      FRACF=DMIN1(FRACF,FARG)
      XMOL(ICAOH2)=FRACF*XWTC(2)    
      IF(XMOL(ICAOH2).GT.XWTC(7)) XMOL(ICAOH2)=XWTC(7)
      XWTC(7)=XWTC(7)-XMOL(ICAOH2)
      XMOL(IFH2O)=(ONE-FRACF)*XWTC(2)
      XMOL(IVH2O)=0.   
   63 CONTINUE
C.....SO FAR, H2O HAS BEEN ASSIGNED TO FREE WATER AND CAOH2, AND
C.....THE CAO INVENTORY HAS BEEN DEFICITED THE MOLAR AMOUNT OF CAOH2    
C.....CREATED.  NOW PARTITION CO2 INVENTORY INTO DOLOMITE AND CACO3.    
      XMIN=DMIN1(XWTC(7),XWTC(8))
      XREQCO=TWO*XMIN 
      IF(XWTC(1).LT.XREQCO) XREQCO=XWTC(1)
      XMIN=HALF*XREQCO
      XMOL(IMCCO3)=XMIN
      XWTC(1)=XWTC(1)-XREQCO    
      XWTC(7)=XWTC(7)-XMIN      
      XWTC(8)=XWTC(8)-XMIN      
C.....PUT BALANCE OF CO2 INTO CALCIUM CARBONATE.    
      XMIN=DMIN1(XWTC(1),XWTC(7))
      XMOL(ICACO3)=XMIN
      XWTC(7)=XWTC(7)-XMIN      
      XWTC(1)=XWTC(1)-XMIN      
C.....ASSIGN REMAINING DATA TO XMOL(I) VECTOR & RENORMALIZE MOLE
C.....FRACTIONS.      
      DO 65 K=6,18    
      IF(K.EQ.7) GO TO 65
      KAS=K-3
      IF(K.GT.6) KAS=K-4
      XMOL(K)=XWTC(KAS)
   65 CONTINUE
      HSUMT=0.D0      
      DO 66 K=1,NMATC 
      HSUMT=HSUMT+XMOL(K)
   66 CONTINUE
      DO 67 K=1,NMATC 
      XMOL(K)=XMOL(K)/HSUMT     
   67 CONTINUE
C.....FIND WEIGHT FRACTIONS H2O, CAOH2, CACO3, AND MGCACO3 IN CONCRETE. 
      HMOLS=ZERO      
      DO 71 K=1,NMATC 
      HMOLS=HMOLS+XMOL(K)*FMMOL(K)
   71 CONTINUE
      XFH2OU=((XMOL(IFH2O)+XMOL(IVH2O))*FMMOL(IFH2O))/HMOLS   
      XFCAOH=(XMOL(ICAOH2)*FMMOL(ICAOH2))/HMOLS     
      XFMGCA=(XMOL(IMCCO3)*FMMOL(IMCCO3))/HMOLS     
      XFCACO=(XMOL(ICACO3)*FMMOL(ICACO3))/HMOLS     
 1134 CONTINUE
C.....DEFINE THE FACTORS EMPLOYED IN THE SHAW VISCOSITY METHODOLOGY     
C.....FOR OXIDE MIXTURES CONTAINING SILICA.
      NMATFI=8
      NMATFF=28
      DO 119 I=1,NMATFF
      XVISC(I)=ZERO   
      SVISC(I)=ZERO   
  119 CONTINUE
      I=IAL2O3
      XVISC(I)=2.0    
      SVISC(I)=6.7    
      I=ITIO2
      XVISC(I)=1.0    
      SVISC(I)=4.5    
      I=IUO2
      XVISC(I)=1.0    
      SVISC(I)=4.5    
      I=IFEO
      XVISC(I)=1.0    
      SVISC(I)=3.4    
      I=IZRO2
      XVISC(I)=1.0    
      SVISC(I)=4.5    
      I=ICR2O3
      XVISC(I)=2.0    
      SVISC(I)=3.4    
      I=IMGO
      XVISC(I)=1.0    
      SVISC(I)=3.4    
      I=ICAO
      XVISC(I)=1.0    
      SVISC(I)=4.5    
      I=IFE2O3
      XVISC(I)=2.0    
      SVISC(I)=3.4    
C.....CALCULATE THE SPECIFIC ENTHALPIES AT THE LIQUIDUS, THE SOLIDUS,   
C.....AND THE TEMPERATURE AT WHICH CONCRETE DECOMPOSITION IS COMPLETE.  
      AL=ZERO
      BL=ZERO
      CL=ZERO
      AS=ZERO
      BS=ZERO
      CS=ZERO
      DNOM=ZERO
      RMASSN=ZERO     
      RMASSD=ZERO     
      DO 199 I=1,NMATC
      XYI=XMOL(I)     
      IF (I .EQ. ICACO3) XYI=ZERO
      IF (I .EQ. ICK2O) XYI=ZERO
      IF (I .EQ. IMCCO3) XYI=ZERO
      IF (I .EQ. ICAOH2) XYI=ZERO
      IF (I .EQ. IFH2O) XYI=ZERO
      IF (I .EQ. ICAO) XYI=XYI+XMOL(ICACO3)+XMOL(IMCCO3)      
     1+XMOL(ICAOH2)   
      IF (I .EQ. IVK2O) XYI=XYI+XMOL(ICK2O)
      IF (I .EQ. IMGO) XYI=XYI+XMOL(IMCCO3)
      IF (I .EQ. IVH2O) XYI=XYI+XMOL(IFH2O)
      AL=AL+XYI*AEQM(I,2)
      BL=BL+XYI*BEQM(I,2)
      CL=CL+XYI*CEQM(I,2)
      AS=AS+XYI*AEQM(I,1)
      BS=BS+XYI*BEQM(I,1)
      CS=CS+XYI*CEQM(I,1)
      DNOM=DNOM+XYI*FMMOL(I)    
      RMASSN=RMASSN+XYI*FMMOL(I)
      RMASSD=RMASSD+XMOL(I)*FMMOL(I)      
  199 CONTINUE
      RDNOM=ONE/DNOM  
      AL=AL*RDNOM     
      BL=BL*RDNOM     
      CL=CL*RDNOM     
      AS=AS*RDNOM     
      BS=BS*RDNOM     
      CS=CS*RDNOM     
      ECL=CL+TCL*(BL+TCL*AL)    
      CCL=BL+TWO*TCL*AL
      ECS=CS+TCS*(BS+TCS*AS)    
      CCS=BS+TWO*TCS*AS
      ECAL=CS+TCAL*(BS+TCAL*AS) 
      RMASS=RMASSN/RMASSD
      RMASSL=RMASS    
C.....CALCULATE THE SPECIFIC ENTHALPIES AT THE TEMPERATURE AT WHICH     
C.....DECOMPOSITION OF CACO3 BEGINS AND THE TEMPERATURE AT WHICH
C.....DECOMPOSITION OF MGCA(CO3)2 IS COMPLETE.      
      A=ZERO
      B=ZERO
      C=ZERO
      DNOM=ZERO
      DO 299 I=1,NMATC
      XYI=XMOL(I)     
      IF (I .EQ. IMCCO3) XYI=ZERO
      IF (I .EQ. ICAOH2) XYI=ZERO
      IF (I .EQ. IFH2O) XYI=ZERO
      IF (I .EQ. ICAO) XYI=XYI+XMOL(IMCCO3)+XMOL(ICAOH2)      
      IF (I .EQ. IMGO) XYI=XYI+XMOL(IMCCO3)
      IF (I .EQ. IVH2O) XYI=XYI+XMOL(IFH2O)
      A=A+XYI*AEQM(I,1)
      B=B+XYI*BEQM(I,1)
      C=C+XYI*CEQM(I,1)
      DNOM=DNOM+XYI*FMMOL(I)    
  299 CONTINUE
      RDNOM=ONE/DNOM  
      A=A*RDNOM
      B=B*RDNOM
      C=C*RDNOM
      ECAS=C+TCAS*(B+TCAS*A)    
      EMCAL=C+TMCAL*(B+TMCAL*A) 
C.....CALCULATE THE SPECIFIC ENTHALPIES AT THE TEMPERATURE AT WHICH     
C.....DECOMPOSITION OF MGCA(CO3)2 BEGINS AND THE TEMPERATURE AT
C.....WHICH DECOMPOSTION OF CA(OH)2 IS COMPLETE.    
      A=ZERO
      B=ZERO
      C=ZERO
      DNOM=ZERO
      DO 399 I=1,NMATC
      XYI=XMOL(I)     
      IF (I .EQ. ICAOH2) XYI=ZERO
      IF (I .EQ. IFH2O) XYI=ZERO
      IF (I .EQ. ICAO) XYI=XYI+XMOL(ICAOH2)
      IF (I .EQ. IVH2O) XYI=XYI+XMOL(IFH2O)
      A=A+XYI*AEQM(I,1)
      B=B+XYI*BEQM(I,1)
      C=C+XYI*CEQM(I,1)
      DNOM=DNOM+XYI*FMMOL(I)    
  399 CONTINUE
      RDNOM=ONE/DNOM  
      A=A*RDNOM
      B=B*RDNOM
      C=C*RDNOM
      EMCAS=C+TMCAS*(B+TMCAS*A) 
      EBWL=C+TBWL*(B+TBWL*A)    
C.....CALCULATE THE SPECIFIC ENTHALPIES AT THE TEMPERATURE AT 
C.....WHICH DECOMPOSTION OF CA(OH)2 BEGINS AND THE RELEASE OF 
C.....FREE WATER IS COMPLETE    
      A=ZERO
      B=ZERO
      C=ZERO
      DNOM=ZERO
      DO 499 I=1,NMATC
      XYI=XMOL(I)     
      IF (I .EQ. IFH2O) XYI=ZERO
      IF (I .EQ. IVH2O) XYI=XYI+XMOL(IFH2O)
      A=A+XYI*AEQM(I,1)
      B=B+XYI*BEQM(I,1)
      C=C+XYI*CEQM(I,1)
      DNOM=DNOM+XYI*FMMOL(I)    
  499 CONTINUE
      RDNOM=ONE/DNOM  
      A=A*RDNOM
      B=B*RDNOM
      C=C*RDNOM
      EBWS=C+TBWS*(B+TBWS*A)    
      EFWL=C+TFWL*(B+TFWL*A)    
C.....CALCULATE THE SPECIFIC ENTHALPY AT THE TEMPERATURE AT WHICH
C.....RELEASE OF FREE WATER BEGINS.
      A=ZERO
      B=ZERO
      C=ZERO
      DNOM=ZERO
      DO 599 I=1,NMATC
      XYI=XMOL(I)     
      A=A+XYI*AEQM(I,1)
      B=B+XYI*BEQM(I,1)
      C=C+XYI*CEQM(I,1)
      DNOM=DNOM+XYI*FMMOL(I)    
  599 CONTINUE
      RDNOM=ONE/DNOM  
      A=A*RDNOM
      B=B*RDNOM
      C=C*RDNOM
      EFWS=C+TFWS*(B+TFWS*A)    
C.....ASSIGN BASEMAT PROPERTIES ACCORDING TO VECTOR NOTATION. 
C.....FIRST SET THE CONCRETE DECOMPSITION TEMPERATURE ACCORDING TO 
C.....THE CONCRETE TYPE....
      IF(ICTC.NE.1) GO TO 2131
      TSCS(1)=1393.0
      GO TO 2135
 2131 CONTINUE
      IF(ICTC.NE.2) GO TO 2132
      TSCS(1)=1403.0
      GO TO 2135
 2132 CONTINUE
      IF(ICTC.NE.3) GO TO 2133
      TSCS(1)=1495.0
      GO TO 2135
 2133 CONTINUE
      TSCS(1)=TDC
 2135 CONTINUE
      TSCS(2)=TSTLS   
      TSCL(1)=TCL     
      TSCL(2)=TSTLL   
      ESCS(1)=ECS     
      ESCS(2)=ESTLS   
      ESCL(1)=ECL     
      ESCL(2)=ESTLL   
      CCSS(1)=CCS     
      CCSS(2)=ESTLS/(TSTLS-298.15D0)      
      CCSL(1)=CCL     
      CCSL(2)=CSTLL   
      CALL DENSS(TSTLS,ROSTLS)  
      CALL DENSS(TSTLL,ROSTLL)  
      DO 88 K=1,16    
      XWTSS(K)=0.D0   
   88 CONTINUE
      XWTSS(2)=1.0    
      IF(NSTEEL.EQ.0.AND.ICTC.LE.4) GO TO 2405 
      IF(NSIMST.EQ.0) GO TO 2405
      DO 2404 K=1,16
      XWTSS(K)=XCSI(K)
 2404 CONTINUE
 2405 CONTINUE
      RETURN
      END