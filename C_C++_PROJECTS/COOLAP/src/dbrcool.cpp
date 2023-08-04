#ifndef DBRCOOL_H
#define DBRCOOL_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <vector>

using namespace std;

#include "meltprop.cpp"
#include "consts.cpp"
#include "coolpropdef.cpp"
#include "wrsttab.cpp"
#include "htpar.cpp"
#include "dhf1d.cpp"
#include "misc_mod.cpp"

/*###############################################################################*/ 
// STRUCT

const int ndrmax = 50, ndzmax=100, ndrndzmax=2500, nthexmax = 10000, nnpsize=50;

struct poolgeom
{
    double hh0=0.0, hhp0=0.0, hhp=0.0, hhf=0.0, rr=0.0, csar=0.0;
};

struct melthex
{
    int nthex = 1;
    double dthex = 0.5;
    double timedq, dqthestep;
    // thex[nthexmax] = {}, qhex[nthexmax] = {}, 
    // vector<double> thex = vector<double>(nthexmax,0.0); // 0:nthexmax
    // vector<double> qhex = vector<double>(nthexmax,0.0); // 1:nthexmax
    // vector<double> qhex = vector<double>(nthexmax,0.0); // 1:nthexmax
    // vector<double> thex = vector<double>(nthexmax+1,0.0); // 0:nthexmax
    double thex[nthexmax+1] = {}; 
    double qhex[nthexmax] = {};
};

struct meltjet
{
    int njet = 50, ibr = 1;
    int jmod = 0; // << juwook 2022
    double dzjet = 1e-2, cbr = 1.0;
    double v0, d0, temp0, vi, di,tempi, vb, db, llj, llbr, zjl, massintot=0.0,eneintot=0.0; 
    double r0, mretot, pin, //! vessel cylinder diameter,total relocated corium mass, pressure difference between in-vessel and containment :: inserted by juwook 220223
            lb=0.0, kk=0.0; //! ICI tube length in apr1400, entrance loss coeff //220227 juk 
    struct mprop mpp;
};
 
struct meltpar
{
    double cvz=0.5, cvr=0.5, ccd=1, chtc=0.5;
    int iranvp=1, itsf=1, id=0,ihtmod;
    int parsdf=1, imod_dmin=1; // << juwook 2022
    // Particle size distributin function (1: R-R 2:T R-R)
    // d_min calculation mode (1: use correlation // 2: use input value) 
    double mass, num, ars, dia, z, r, vz, vr,dmmb, tsf;
    double deltd=0.0, dmin=0.0, qfluxmean=0.0; // << juwook 2022
    // minimum diameter of particle - 0414 juwook //mean flux
    struct mprop mpp;
};

struct dbstack
{
    int ndz = 0; 
    vector<int> ihtmod = vector<int>(ndzmax,0);
    double zb[ndzmax+1], zc[ndzmax], dz[ndzmax], // ! z node boundary positions, sizes
    //   por[ndzmax],  // ! local porosity [0.5 by default, can be changed by input or model]
      mass[ndzmax], num[ndzmax], ars[ndzmax], // ! mass,particle number,surface area per node 
      qrel[ndzmax], tsf[ndzmax], dia[ndzmax]; // ! heat release per node [W], dia., surface temp. 
    vector<double> por = vector<double>(ndzmax, 0.0);
    struct mprop mpp[ndzmax];
};

struct debrbed // << by juwook 2022
{
    int ndr=20,  // ! number of nodes (<ndrmax)
        idiv=2,  // ! division scheme: 1: even floor area, 2: even distance (default)
        irep=1,  // ! 1:consider avalanche with fixed repose angle (default), 0:not
        irepinlate=1, // ! consider late in-side avalanche 0:no(default), 1:yes
        idhf=1,     // ! 1:use Rahman(2013) 1D top flooding DHF model as heat transfer limit (default),
                    // ! 2:use Reed(1982) model instead of Rahman,
                    // ! 0:not use DHF model (particle heat transfer based eval.)

        ishape=1,   // Debris bed shape model option
                    // 1:Original model
                    // 2:Eunho Kim's model
        maxalphanode=0;
    double  dzdb=5e-2,  // ! standard z node size (m)
            angrep=10,  // ! avalanche repose angle (deg)(default 10 deg, relatively flat bed
                        // ! considering boiling agitation), 
            tangrep,    // ! tangrep=tan(angrep) for comparison
            angrepin=35, tangrepin, // ! angrep for late in-side avalanche (for irepinlate=1)
            chtc=0.1, // ! heat transfer factor for particles in debris bed
            epor=0.5, // ! porosity, default 0.5 based on FARO exp. observations 
            fdhf=1,   // ! factor for DHF (to consider 2D enhance effect)
            qrel,qrelbtm,qdc,   // ! heat release and decay heat power (W)
            qreltot=0, qdctot=0, // ! total released heat and decay heat (J)   
            masstolmp=0, // ! total mass moved to lump (remelt) (kg)
            alpinidhf=0, // ! user initial guess void frac. for Rahman's DHF model (non-zero => use it)
            vlbtm=0,    // ! water inlet superficial velocity (m/s) at bottom

            Vtot=0.0,     // total volume of debris bed (m^3)
            vp=0.0,  // particle falling velocity (m/s)//conical
            tau=1.0, 
            massseive[ndrmax], 
            dbrepose=0.0, dbheight=0.0, dbradius=0.0, 
            jj_g=0.0, mdot, cdpc, 
            x1,x2,x3,x4,x5,x6,x7,x8;
    double rb[ndrmax+1],rc[ndrmax], dar[ndrmax]; // ! r node boundary/center positions, bottom area
    struct dbstack dbr[ndrmax];
};

struct meltlmp
{
    int itcr=1; // ! particle temperature criterion for merge to lump
      // ! 0: no merge
      // ! 1(default): average temperature criterion tav>tmel+ftcr*0.5*(tliq-tsol), ftcr=-1..1
      // ! 2(default): surface temperature criterion tsf>tsol+ftcr*(tliq-tsol), ftcr=-1..1
    double ftcr = 1,  // ! factor for the merge criterion (see above)
         fhtlmp = 1, fqbdhf=0.5; // ! factor for consideration of melt lump top heat release
         // ! q"btm = fhtlmp*lamm*(Tm-Tsat)/(Hlmp/2), q"btm = min(q"btm, fqbdhf*DHF)
    double mass, thk, qreltot=0, qdctot=0, qrel,qdc;
      // ! mass (kg), thickness for heat transfer eval. (m), 
      // ! total released heat and decay heat (J), heat release and decay heat power (W)
    struct mprop mpp; // ! melt lump physical properties
};

struct parsieve
{
    int nn=19; // ! element no. of diab
    double diab[19] = {5e-5, 7.5e-5, 1e-4, 1.5e-4, 2.3e-4, 3.4e-4, 5e-4, 7.5e-4, 1e-3,
        1.5e-3, 2.3e-3, 3.4e-3, 5e-3, 7.5e-3, 1e-2, 1.5e-2, 2.3e-2, 3.4e-2, 5e-2}, 
       mass[20] = {0}, masstot = 0.0;
};

struct meltdebr
{
    struct meltjet mj; // ! melt jet
    struct meltpar mp; // ! melt particle
    struct debrbed db; // ! debris bed
    struct meltlmp ml; // ! melt lump
    struct melthex hx; // ! particle heat change buffer
    struct timetab vjin, djin, tjin; // ! melt inlet condition
    struct mpropbase mppb; // ! melt physical properties base data
                            // ! (melt properties are included in each component)

    struct poolgeom plg; // ! pool geometry
    struct coolprop clp; // ! coolant physical properties
    int ic; // ! index of "can" to which mbrq model is conected
                  // ! if 0 is given, mbrq is not implemented.
    struct parsieve psievegen, psieve;
                  // ! particle sieve record for fresh (just generated) particles
                  // ! and for solid sediment particles
    // ! parameters for particle size generater with default numbers given here
    int npsize=20,  // ! actually used sample size
      ipsize=1,    // ! particle size option: 1=size distribution correlation (default),
                    // !   2=distribution with given dmm, 3:fixed size
      ipsizeran=0;   // ! option 0: not use random sampling but use fixed sizes 
                    // ! at center of layers (default), 1: use random sampling (LHS)
    double cdmm=1.0,  // ! modification factor for mass median diameter used if ipsize=1
      cpsize=3e-3; // ! given dmm if ipsize=2 or fixed particle size if ipsize=3

    // ! output and output time control variables
    int ujh=21,ulh=22,udh=23,uhh=24,ubh=25,udl=26,upr=27,upt=28,ups=29;
         // ! output units for jet history, lump history, debris bed history, hex history,
         // ! balance history, debris bed list, particle record, particle tracking,
         // ! particle sieve, with their default numbers
    double dtouth,dtoutl,dtoutp, touth=0.0,toutl=0.0,toutp=0.0;
         // ! time steps for output of histories, lists and particle, init. vars for output times
         // ! setting dt*** = 0 supresses the corresponding output
    int noutp; // ! step number to end the particle data output
  
  // ! some common variable for mass/energy balance check
    double masspartot = 0.0, // ! a buffer to catch total particle mass (kg)
      wpoolmass, wpoolmassevap,wpoolenth,wpooleneevap, // ! water pool mass, energy totals (J)
      wpoolmass0, enthvap0, // ! initial value of water mass, vapor enthalpy
      totqsrcex=0.0; // ! total external heat source read by hotcan (J) (mbrq+body)
};

// ! decay heat related parameters
double t0sdown=0.0, powopr=0.0, coremass=0.0, fqdcex=0.0;
    // ! time after shutdown at time=0 (s), operation power (W), core debris mass (kg),
    // ! fraction of decay heat outside the coremass (e.g. volatilized before melt production)
    // ! these parameters are confined inside module, setting is done by initdecayheat(...)
double qdcfrac=0.0, qdecay=0.0, qdecaymass=0.0;
    // ! decay heat rate by fraction to operation power, 
    // ! decay heat power (W), decay heat power per mass of core debris (W/kg)
    // ! these are calculated by decayheat(...) for given time
int iresdh=1; // ! include residual decay heat (the part possessed by mass not
    // ! in the melt/debris model, e.g. remaining in reactor vessel) 1: include, 0: not
    // ! (decay heat in the jet/particles is still neglected, it is only a small portion)
/*###############################################################################*/ 








/*###############################################################################*/ 
// RANDOM FUNCTION :: WHANG
double randomX()
{
    double x;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0,1);
    x = dis(gen);
    return x;
}
/*###############################################################################*/ 
 






















/*###############################################################################*/ 
void depomlumpjet(meltjet mj, meltlmp &ml, mpropbase mppb, double dt)
// ! deposition of direct melt jet to melt lump
{
    double ma, m, e; 

    ma = mj.vb * pi/4.0*pow (mj.db,2) * mj.mpp.rhom * dt; //  ! mass addition
    m = ml.mass + ma;   
    e = (ml.mpp.enem*ml.mass + ma * mj.mpp.enem)/m;
    ml.mass = m;
    ml.mpp.enem = e;

    menetemp(mppb, ml.mpp.enem, ml.mpp.tm);
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void parsizegen (meltdebr &md, double &dmmsv)
{
    static int istpsz=0,  // ! status of particle size generator 1:initialized, 0:not
        countpsz;             // ! counter to control the rolling (stores the index of datum previously used)
    static double ffpsz[nnpsize+1],  // ! sections of cumulative probability 0--1
        ddpsz[nnpsize];    // ! particle diameters for the sampled fractions
    double frpsz[nnpsize];  // ! sampled ramdom number [0..1] [if no LHS, always 0.5]

    int i, j;
    double x, drl, drt, tsub, nsub, rrrho, nnrho, fsb, dmm;
    static double dffpsz;
    double qfluxmean,a=33.7,b=0.0408, hlatent=0.0,x1,x2; // << juwook 2022
    // save :: ffpsz, dffpsz, ddpsz, istpsz, countpsz

    // ! rayleigh instability droplet size as max. limit droplet size
    drl = md.mj.di*1.89;
    // if(drayleigh !=-1){drayleigh = drl;}
    // ! bond number based (R-T instability) max. limit droplet size (empirical)
    drt = 10.0*sqrt(md.mj.mpp.sigm/grav/md.mj.mpp.rhom);
    // if(draytay !=-1){draytay = drt;}

    // ! just use given fixed size
    if(md.ipsize == 3)
    {
        md.mp.dia = md.cpsize;
        // if(dmmsv !=-1){dmmsv = md.mp.dia;} // ! preserve evaluated dmm for output
        dmmsv = md.mp.dia;
        return;
    }
    // end if

    // ! initialization of particle size bins
    if(istpsz == 0){
        ffpsz[0] = 0.001;       // ! sample in 0.1%-99.9% CDF
        ffpsz[md.npsize] = 0.999;
        dffpsz = (ffpsz[md.npsize]-ffpsz[0])/(md.npsize);
        for (i=1; i <= md.npsize-1; i++){
            ffpsz[i] = ffpsz[i-1]+dffpsz;}
        // end do
        istpsz = 1; // ! set flag initialized.
        countpsz = md.npsize; // ! set counter to be the end val. to force data gen. at 1st call
    }
    // end if

    if(countpsz < md.npsize)
    {
        // ! when size data is available, increment the magazine and give the next datum
        countpsz = countpsz + 1;
        md.mp.dia = ddpsz[countpsz-1];
    }
    else
    {
        // ! generate a new sample of the sizes
        // ! mass median diameter from exp. correlation
        if (md.ipsize == 1) // ! dmm correlation
        {
            tsub = max((md.clp.tspt - md.clp.tl), 0.0);
            nsub = md.clp.cpl*tsub/md.clp.hfg;
            rrrho = md.clp.rhov/md.clp.rhol;
            nnrho = md.mj.mpp.rhom/md.clp.rhol;
            fsb = 25*(2.2-exp(-13*nsub));
            dmm = fsb*pow(rrrho,0.3333)/pow(nnrho,0.6667)*sqrt(md.mj.mpp.sigm/grav/md.mj.mpp.rhom);
            dmm = md.cdmm*dmm; // ! put an arbitrary factor for parametric tests
        }
        else // ! given dmm (ipsize == 2)
        {
            dmm = md.cpsize;
        }
        // endif
        // if(dmmsv!=-1) {dmmsv = dmm;} // ! preserve evaluated dmm for output
        dmmsv = dmm;
    
        // ! sample ffpsz() and get ddpsz()=particle size sample
        for (i=1; i<=md.npsize; i++) // do i = 1, md.npsize
        {
            if(md.ipsizeran == 0) {x = 0.5;} // ! just take the center of the F section
            else
            {
                x = randomX(); // ! stratified random sample (LHS)
            }
                // end if
                frpsz[i-1] = x*dffpsz + ffpsz[i-1];
                // << juwook 2022
                if (md.mp.parsdf==1) {
                    ddpsz[i-1] = dmm*pow((-log(1.0-frpsz[i-1])/log(2.0)),(0.6667)); // ! droplet dia. by Rosin-Rammler distribution with n=1.5
                }
                else if (md.mp.parsdf==2) {
                    if (md.mp.imod_dmin==1) {
                        qfluxmean=md.mp.qfluxmean;
                        hlatent=max(md.clp.hfg, 0.0);
                        x1=pow((pow(qfluxmean,(4.0))+8.0*a*b*md.clp.rhov*(md.mj.mpp.rhom-md.clp.rhol)*grav*pow(hlatent,(3))*md.clp.muv*qfluxmean),0.5);
                        x2=2.0*b*md.clp.rhov*(md.mj.mpp.rhom-md.clp.rhol)*grav*pow(hlatent,(2));
                        md.mp.dmin=max((pow(qfluxmean,(2.0))+x1)/x2,0.0);
                        md.mp.dmin=min(md.mp.dmin,dmm);
                        ddpsz[i-1] = pow((pow(dmm,(1.5))-pow(md.mp.dmin,(1.5))),(0.6667))*pow((-log(1.0-frpsz[i-1])/log(2.0)),(0.6667))+md.mp.dmin;
                    }
                    else if (md.mp.imod_dmin==2) {
                        md.mp.dmin=min(md.mp.dmin,dmm);
                        ddpsz[i-1] = pow((pow(dmm,(1.5))-pow(md.mp.dmin,(1.5))),(0.6667))*pow((-log(1.0-frpsz[i-1])/log(2.0)),(0.6667))+md.mp.dmin; // !droplet dia. by Truncated R-R distribution with n=1.5
                    }
                    // endif
                }
                // endif
            ddpsz[i-1] = min(min(ddpsz[i-1], drl), drt);
        }
        // end do
    
        // ! shaffle the sample 
        for (i= 1; i <= md.npsize; i++) // do i = 1, md.npsize
        {
            x = randomX();
            frpsz[i-1] = x;
        }
        // end do
        for (i=1; i<md.npsize; i++) // do i = 1, md.npsize-1
        {
            for (j=i+1; j<=md.npsize;j++) // do j = i+1, md.npsize
            {
                if(frpsz[i-1] > frpsz[j-1]){
                    x = frpsz[i-1];
                    frpsz[i-1] = frpsz[j-1];
                    frpsz[j-1] = x;
                    x = ddpsz[i-1];
                    ddpsz[i-1] = ddpsz[j-1];
                    ddpsz[j-1] = x;
                }
                // end if
            }
            // end do
        }
        // end do
    
        // ! give out the 1st datum
        countpsz = 1;
        md.mp.dia = ddpsz[countpsz-1];
    }
    // end if
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void recparsieve(meltpar &mp, parsieve &psv) 
{
    double mass, dia;
    int i;

    dia = mp.dia;
    mass = mp.mass*mp.num;
    psv.masstot = psv.masstot + mass;
    for (i=1; i<=psv.nn; i++) // do i = 1, psv.nn
    {
        if(dia <= psv.diab[i-1])
        {
            psv.mass[i-1]=psv.mass[i-1]+mass;
            return;
        }
        // end if
    }
    // end do
    psv.mass[psv.nn]= psv.mass[psv.nn] + mass;
    // return
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void depomhex(melthex &mh, double timenow, double timesto, double qin)
{
    int i, inow;
    double qrate; 

    if (timesto <0.0 || timesto > mh.thex[nthexmax]){
        printf ("'depomhex{ out of thex range.'");
        exit(EXIT_FAILURE);
    }
    i = ceil(timesto/mh.dthex);

    mh.nthex = max(i, mh.nthex);
    qrate = qin/mh.dthex;
    mh.qhex[i-1] = mh.qhex[i-1] + qrate;

    if(timenow > mh.timedq){
        mh.dqthestep = 0.0;
        mh.timedq = timenow;
    }
    // end if

    // ! during a time step (for melt and coolant evolution), the increased heat rate 
    // ! for the time span "begining of the time step to now" is igored when qhex is
    // ! read. the part of increase is stored for the usage for compensation.
    // ! this buffer "dqthestep" should be cleared when a new time step begins.
    // ! it is done above or in "readmhex".
    inow = ceil(timenow/mh.dthex);
    if(i == inow){
        mh.dqthestep = mh.dqthestep + qrate; // ! accumulate heat rate addition in the time step
    }
    return; 
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void depomlump(meltpar mp,meltlmp &ml,mpropbase mppb,double time)
{
    /* 
    ! deposit melt particles to melt lump
    ml: inout 
    */
    double ma, m, e;

    ma = mp.mass*mp.num;
    m = ml.mass + ma;
    e = (ml.mpp.enem*ml.mass + ma*mp.mpp.enem)/m;
    ml.mass = m;
    ml.mpp.enem = e;
    menetemp(mppb, ml.mpp.enem, ml.mpp.tm);

    // # ! no other physical properties are calculated.
    // # ! they are needed when further modeling will be done.
    return;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void depodbed(meltpar mp,debrbed &db,mpropbase mppb,double time)
{
    /* 
    ! deposit melt particles to debris bed
    db:: inout 
    */
    int i, iz, ia;
    int iseive;
    double ma, m,e,n, tang, va,da,sa,s0,v0;

    // ! find the radial position
    if(mp.r < 0.0 || mp.r > db.rb[db.ndr]){
        /* """LATER""" */
        printf("'call errorhdr'"); //call errorhdr(mess="depodbed: out of radius range.",act="stop")
        exit(EXIT_FAILURE);
    }
    // end if
    for (i=1; i<=db.ndr; i++){ // do i = 1, db.ndr
        if(mp.r <= db.rb[i]){break;}
    }
    // end do
    ia = i; // ! node where a particle is arriving
    // by Juwook 2022
    iseive = i; // node num for seive

    // ! check avalanche at angle > db.angrep (comparison by tangent)
    /* """LATER: CHECK THE RANGE, dbr(i)""" */
    if(db.irep > 0){  // ! ia = index of the arriving node
        for (i=ia; i>=2; i--){ // // do i = ia, 2, -1 // ! check fall to inside
            tang = (db.dbr[i-1].zb[db.dbr[i-1].ndz]-db.dbr[i-2].zb[db.dbr[i-2].ndz])
                /(db.rc[i-1]-db.rc[i-2]);
            if(tang < db.tangrep){break;}
        }
        // end do
        
        if(i < 2){ia = 1;}
        else{ia = i;}
        // end if

        for (int ii=ia; i<=db.ndr; i++) { // do i = ia, db.ndr-1 // ! check fall to outside
            tang = (db.dbr[i-1].zb[db.dbr[i-1].ndz]-db.dbr[i].zb[db.dbr[i].ndz]) 
                /(db.rc[i]-db.rc[i-1]);
            if(tang < db.tangrep){break;}
        }
        // end do

        if(i > db.ndr-1){ia = db.ndr;}
        else{ia = i;}
        // end if
    }
    // end if

    i = ia;
    // ! deposit the particle to debris stack at i
    if(db.dbr[i-1].dz[db.dbr[i-1].ndz-1] > db.dzdb){ // ! if the top node exceeded standard size,
        if(db.dbr[i-1].ndz < ndzmax){                     // ! add a vertical node and initialize it
            db.dbr[i-1].ndz = db.dbr[i-1].ndz + 1;
            iz = db.dbr[i-1].ndz;
            db.dbr[i-1].dz[iz-1] = 0.0;
            db.dbr[i-1].mass[iz-1] = 0.0;
            db.dbr[i-1].ars[iz-1] = 0.0;
            db.dbr[i-1].num[iz-1] = 0.0;
            db.dbr[i-1].tsf[iz-1] = 0.0;
            db.dbr[i-1].dia[iz-1] = 0.0;
            db.dbr[i-1].qrel[iz-1] = 0.0;
            db.dbr[i-1].ihtmod[iz-1] = 0;
            db.dbr[i-1].mpp[iz-1].enem = 0.0;
        }
        else{
            iz = db.dbr[i-1].ndz;
            /* """LATER""" */
            printf("'call errorhdr\n'"); // call errorhdr(mess="depodbed: ndz >= ndzmax // !",act="cont")
            exit(EXIT_FAILURE);
        }
        // end if
    }
    else{iz=db.dbr[i-1].ndz;}
    // end if

    // ! mass, internal energy, number conservation
    // ! diameter and surface area can change by temperature change
    ma = mp.mass*mp.num;
    m = db.dbr[i-1].mass[iz-1] + ma;
    e = (db.dbr[i-1].mpp[iz-1].enem*db.dbr[i-1].mass[iz-1] + ma * mp.mpp.enem)/m;
    // ! new temp and phys. props.
    menetemp(mppb,e,db.dbr[i-1].mpp[iz-1].tm);
    getmprop(mppb,db.dbr[i-1].mpp[iz-1]);

    // by Juwook 2022
    db.massseive[i-1] = db.massseive[i-1]+ma;


    // ! volume and surface area at new temp.
    // ! receiving node
    v0 = db.dbr[i-1].mass[iz-1]/db.dbr[i-1].mpp[iz-1].rhom;
    s0 = pow((36.0*pi*db.dbr[i-1].num[iz-1]*pow(v0,2)),(1.0/3.0));
    // ! arriving particles
    va = ma/db.dbr[i-1].mpp[iz-1].rhom;
    sa = pow((36.0*pi*mp.num*pow(va,2)),(1.0/3.0));

    // ! mass, energy update for the receiving node
    db.dbr[i-1].mass[iz-1] = m;
    db.dbr[i-1].mpp[iz-1].enem = e;

    // ! conserved surface area, sauter mean dia. at new temperature
    db.dbr[i-1].ars[iz-1] = s0 + sa;
    db.dbr[i-1].num[iz-1] = pow(db.dbr[i-1].ars[iz-1],3)/pow((v0+va),2)/36.0/pi;
    db.dbr[i-1].dia[iz-1] = 6.0*(v0+va)/db.dbr[i-1].ars[iz-1];

    // ! surface temperature: surface area weighted average
    db.dbr[i-1].tsf[iz-1] = (db.dbr[i-1].tsf[iz-1]*s0 + mp.tsf*sa)/db.dbr[i-1].ars[iz-1];
    db.dbr[i-1].tsf[iz-1] = min(db.dbr[i-1].tsf[iz-1], db.dbr[i-1].mpp[iz-1].tm);

    // ! node size, dz
    db.dbr[i-1].dz[iz-1] = (v0+va)/db.dbr[i-1].por[iz-1]/db.dar[i-1];

    // ! porosity is, right now, kept as default value, 0.5; can be changed at initialization

    return;
}
/*###############################################################################*/ 





































/*###############################################################################*/ 
void trackpar (meltdebr &md, double time, int iout, double &kepar)
{
    double v0, v1, t, dt, u, w, r, z, m, nstep=50.0, 
        // ! nstep: standard step number to calculate till settlement
        vr, re, cd, ff, fz, fr, tav,tsf, rhovf,cpvf,lamvf,muvf,am, qflx, hrel, hreltot, 
        x0,x1,x2,x3,x4, dbndmx, dbnd,htc,t0, tvf, zgen,tage;
    int i, i1, imrg;

    // whang add for opt declare
    static bool initCall = true;

    ofstream opt, opr;
    opt.precision(4);
    opr.precision(4);
    
    if (initCall){
        opt.open("mbrq/opt", ios::out);
            opt << "" << endl;
        opt.close();

        opr.open("mbrq/opr", ios::out);
            opr << "# particle summary records" << endl;
            opr << "#1:id 2:tgen 3:zgen 4:tset 5:rset 6:irset 7:mass 8:num 9:diaset 10:eneset 11:tempset 12:tsfset" << endl;
        opr.close();

        initCall = false;
    }
    // whang add for opt declare



    recparsieve(md.mp, md.psievegen); // ! sieve fresh particles and rec.

    t = time;
    r = md.mp.r;  
    z = md.mp.z;
    u = md.mp.vr; 
    w = md.mp.vz; // + random
    m = md.mp.mass;

    // ! rough guess for terminal velocity (newton refime, cd=0.44 giving fastest one)
    v0 = sqrt(4*md.mp.dia*grav*(md.mp.mpp.rhom-md.clp.rhol)/(3*md.clp.rhol*0.44));
    v1 = fabs(w); // ! initial velocity
    dt = min( z/min(v0,v1)/nstep, 0.1);

    // ! inital value of tsf=tav=melt jet temp.
    tav = md.mp.mpp.tm;
    tsf = tav;
    md.mp.tsf = tsf;   // ! initialize tsf for the 1st // write
    md.mp.ihtmod = 0;  // ! initialize heat transfer mode buffer for the 1st // write

    // ! particle data output
    // TMP
    if(iout > 0){
        opt.open("mbrq/opt",ios::app);
        opt << "#id= " << md.mp.id << endl;
        opt << "#1:time 2:tage 3:z 4:r 5:vz 6:vr 7:mass 8:num 9:dia 10:ene 11:temp 12:tsf 13:ihtmod" << endl;
    }
    // end if
    tage = 0; 
    zgen = z; // ! reset particle age, set generation height

    //TMP
    // if(time>=0.4){
    //     printf("WT & z = %f\n", z);
    // }


    while(z>0) // do
    {
        // ! particle tracking data output
        if(iout > 0){
            opt << "\t" << fixed << time << "\t" << tage << "\t" << scientific << z << "\t" << r << "\t" << w << "\t" 
                << u << "\t" << m << "\t" << md.mp.num << "\t" << md.mp.dia << "\t" << md.mp.mpp.enem << "\t" 
                << md.mp.mpp.tm << "\t" << md.mp.tsf << "\t" << fixed << md.mp.ihtmod << endl;
        }


        // ! particle kinetics
        vr = sqrt(pow(u,2) + pow(w,2));
        re = max(md.clp.rhol*vr*md.mp.dia/md.clp.mul, 1.0);
        cd = max(max(24/re, 18.5/pow(re,0.6)), 0.44);
        ff = md.mp.ccd * cd * 0.25*pi*pow(md.mp.dia,2) * 0.5*md.clp.rhol*vr;
        // !    fz = -ff*w - grav*(mp.mpp.rhom-clp.rhol)*pi/6*mp.dia**3
        // !    fr = -ff*u
        // !    u = u + dt*fr/m
        // !    w = w + dt*fz/m
        u = u/(1 + dt*ff/m); // ! more stable by implicit scheme on velocity
        w = (w - dt/m*grav*(md.mp.mpp.rhom-md.clp.rhol)*pi/6*pow(md.mp.dia,3))/(1 + dt*ff/m);
        r = r + dt*u;
        z = z + dt*w;

        // ! check if radial position exceeds
        if(r < 0){
            r = -r;
            if(u < 0) {u = -u;}
        }
        else if (r>md.plg.rr){
            r = 2*md.plg.rr - r;
            if ( u > 0) {u = -u;}
        }
        // end if

        // ! vapor film properties need evaluation here (at ptot,0.5*(tspt+tsf))
        tvf = max(0.5*(tsf+md.clp.tspt),md.clp.tspt);
        wrsteamtab(md.clp.ptot,tvf,1, rhovf,x0,x1,x3,x4,cpvf,i1);
        lamvf = tconwaterf(rhovf,tvf);
        muvf = viscwaterf(rhovf,tvf);

        // ! heat transfer
        htranspar(tsf,md.mp.dia,md.mp.mpp.epsm, 0,vr,md.clp.ptot,0, 
            md.clp.tspt,md.clp.tl,md.clp.tspt,md.clp.hfg,md.clp.rhov,md.clp.rhol,md.clp.cpv,md.clp.cpl,
            md.clp.muv,md.clp.mul,md.clp.lamv,md.clp.laml,md.clp.sig, 
            rhovf,cpvf,lamvf,muvf, 1, 
            i1, qflx, x0, x1, x2);
        qflx = qflx * md.mp.chtc; // ! modify heat flux by given factor
        md.mp.ihtmod = i1;
        // !   hrel = qflx*mp.ars*dt // ! heat release per particle

        // ! energy update by an implicit temp scheme for stability
        x0 = fabs(qflx/max(tav-md.clp.tl,10.0)); // ! heat transfer coeff.
        x1 = dt*md.mp.ars*x0;
        x2 = x1/md.mp.mass/md.mp.mpp.dedt;
        hrel = x1*(tav - md.clp.tl)/(1+x2); // ! heat release per a particle

        // ! heat exchange deposit
        hreltot = hrel * md.mp.num;  // ! total heat release
        depomhex(md.hx, time, time+tage, hreltot); // ! deposit in heat exchange buffer

        // ! energy conservation per particle
        md.mp.mpp.enem = md.mp.mpp.enem - hrel/md.mp.mass; // ! energy reduce per mass
        md.mp.mpp.enem = md.mp.mpp.enem + qdecaymass*dt; //  decay generation in particle -0414 by Juwook 2022

        menetemp(md.mppb,md.mp.mpp.enem,md.mp.mpp.tm); // ! update temperature
        tav = md.mp.mpp.tm;

        // ! particle properties at updated temperature
        getmprop(md.mppb,md.mp.mpp);
        am = md.mp.mpp.lamm/md.mp.mpp.cpm/md.mp.mpp.rhom;
        md.mp.dia = pow((6/pi*md.mp.mass/md.mp.mpp.rhom),(1.0/3.0));
        md.mp.ars = pi*pow(md.mp.dia,2);

        // ! consideration of surface temperature for heat transfer
        // ! am : thermal diffusivity of melt
        // ! dbnd: thermal boundarylayer thickness is the particle
        // ! dbndmx=dia/2 : max guess
        if(md.mp.itsf == 0) // ! just use tav for heat transfer
        {
            tsf = tav;
        }
        else    // ! consider surface temperature drop of particles (itsf==1, default)
        {
            dbndmx = max(md.mp.dia/2,1e-4); // ! max limit of boundary layer = particle radius
            if(dbnd <= 0){dbnd = sqrt(6*am*dt);}
            else
            {
                if(dbnd < dbndmx)
                {
                x0 = dbnd/dbndmx;
                x1 = 3*am*dt/max(dbnd*(1-0.75*x0+0.2*pow(x0,2)),1e-6);
                x1 = max(x1, 0.0);
                dbnd = dbnd+x1;
                } // end if
                dbnd = min(dbnd, dbndmx);
            }
            // end if
            x0 = dbnd/dbndmx;
            t0 = md.clp.tspt;
            htc = fabs( qflx/max(tsf-t0, 1.0) );
            x1 = max(1-x0*(1-0.5*x0+0.1*pow(x0,2)),0.0);
            x2 = max(dbnd*htc/md.mp.mpp.lamm/2,0.0);
            tsf = (tav+x1*x2*t0)/max(1+x1*x2,1e-6);
            tsf = min(max(tsf,(t0+tav)/2),tav);
        }
        // end if
        md.mp.tsf = tsf; // ! save calculated surface temp.
        tage = tage + dt;
        // by Juwook 2022
        if (z <= 0.0){
            kepar = pow(u,2) + pow(w,2);
            break;
        }
        // if (z <= 0) {break;}
    } // end do
    // while(z>0);

    // ! save the settling position
    md.mp.z = z;
    md.mp.r = r;

    // ! particle data output
    if(iout > 0)
    {
        // opt.open("mbrq/opt", ios::app);
        opt << "\t" << fixed << time << "\t" << tage << "\t" << scientific << z << "\t" << r << "\t" << w << "\t" 
            << u << "\t" << m << "\t" << md.mp.num << "\t" << md.mp.dia << "\t" << md.mp.mpp.enem << "\t" 
            << md.mp.mpp.tm << "\t" << md.mp.tsf << "\t" << fixed << md.mp.ihtmod << endl;
        opt << "" << endl;
        opt.close();

        for (int i = 1; i <= md.db.ndr; i++) // do i = 1,md.db.ndr
        {
            if(md.db.rb[i] >= r) {break;}
        } // end do
        opr.open("mbrq/opr", ios::app);
        opr << "\t" << fixed << md.mp.id << "\t" << time << "\t" << scientific << zgen << "\t" << fixed << tage << "\t" << scientific << r << "\t" 
            << fixed << i << "\t" << scientific << md.mp.mass << "\t" << md.mp.num << "\t" << md.mp.dia << "\t" << md.mp.mpp.enem << "\t" 
            << md.mp.mpp.tm << "\t" << md.mp.tsf << endl; 
        opr.close();
    }
    // end if

    // ! settlement
    // ! check merge criterion and deposite to lump or debris bed
    imrg = 0; // ! flag for merge
    x0 = md.ml.ftcr*(md.mp.mpp.tliq-md.mp.mpp.tsol); // ! ftcr=-1..1
    if (md.ml.itcr == 2) // ! surface temp. criterion
        {if(tsf > md.mp.mpp.tsol + x0) {imrg = 1;}}
    else if(md.ml.itcr == 1)// ! average temp. criterion (default)
        {if(tav > md.mp.mpp.tmel + 0.5*x0) {imrg = 1;}}
    else {imrg = 0;} // ! no merge (itcr == 0)
        
    // end if
    if(imrg == 1){
        depomlump(md.mp,md.ml,md.mppb,time);
    }
    else{
        // printf("trackpar-tmp2\n");
        recparsieve(md.mp, md.psieve); // ! sieve solid particles and rec.
        depodbed(md.mp,md.db,md.mppb,time);
    }
    // printf("trackpar-tmp3\n");
    // end if
    return; 
}
/*###############################################################################*/ 











































/*###############################################################################*/ 
void initdbed(debrbed &db, double rr)
// ! initialize debris bed
{
    double e, drev;
    int i; 

    // ! descritize the area for rings of uniform areas and initialize vertical stacks
    db.rb[0] = 0.0;
    drev = rr/db.ndr;
    for (i=1; i<=db.ndr; i++) // do i = 1,db.ndr
    {
        if(db.idiv == 1) { // ! even floor area per node
            db.rb[i] = sqrt(pow(db.rb[i-1],2)+pow(rr,2)/db.ndr); // ! node boundary position, rb
        }
        else { // ! idiv == 2, even distance
            db.rb[i] = db.rb[i-1] + drev; // ! node boundary position, rb
        }
        db.rc[i-1] = 0.5*(db.rb[i-1] + db.rb[i]);    // ! node center position, rc
        db.dar[i-1] = pi*(pow(db.rb[i],2) - pow(db.rb[i-1],2)); // ! bottom area
        db.dbr[i-1].ndz = 1; // ! clear vertical debris stack
        db.dbr[i-1].zb[0] = 0.0; // ! bottom position
        db.dbr[i-1].dz[0] = 0.0 ;// ! empty 1st layer
        db.dbr[i-1].mass[0] = 0.0;
        db.dbr[i-1].ars[0] = 0.0;
        db.dbr[i-1].num[0] = 0.0;
        db.dbr[i-1].tsf[0] = 0.0;
        db.dbr[i-1].dia[0] = 0.0;
        db.dbr[i-1].qrel[0] = 0.0;
        db.dbr[i-1].mpp[0].enem = 0.0;
        fill(db.dbr[i-1].ihtmod.begin(), db.dbr[i-1].ihtmod.end(),0);
        fill(db.dbr[i-1].por.begin(),db.dbr[i-1].por.end(),db.epor);
    }

    // by Juwook 2022
    // ! repose angle tangent as comparison index
    // ! auto detect the data in radian and conver to degree
    // if(db.angrep < 0.80) { // ! it is radian (0.8rad = 45.8deg)
    //     db.angrep = db.angrep/pi*180.0;
    // }
    db.tangrep = tan(db.angrep*pi/180.0);
    // ! for late in-side avalanche
    // if(db.angrepin < 0.80) { // ! it is radian (0.8rad = 45.8deg)
    //     db.angrepin = db.angrepin/pi*180.0;
    // }
    db.tangrepin = tan(db.angrepin*pi/180.0);
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void fixdbednode(debrbed &db)
{
    // ! fix node bounds (zb,zc) considering updated node contents (vertical sizes)
    int i, iz;

    for (i=1; i<=db.ndr; i++){ // do i = 1, db.ndr
        for (iz=1; iz<=db.dbr[i-1].ndz; iz++){ // do iz = 1, db.dbr[i].ndz
            if(db.dbr[i-1].mass[iz-1] > 0.0){
                db.dbr[i-1].dz[iz-1] = db.dbr[i-1].mass[iz-1]/db.dbr[i-1].mpp[iz-1].rhom/db.dbr[i-1].por[iz-1]/db.dar[i-1];
            }
            else{
                db.dbr[i-1].dz[iz-1] = 0.0;
            }
            // end if
            db.dbr[i-1].zb[iz] = db.dbr[i-1].zb[iz-1] + db.dbr[i-1].dz[iz-1];
            db.dbr[i-1].zc[iz-1] = db.dbr[i-1].zb[iz-1] + db.dbr[i-1].dz[iz-1]/2.0;
        }
        // end do
    }
    // end do

    return; 
}
/*###############################################################################*/ 



void conical(meltdebr &md,double dt) // !updated 0523_juwook
{
  double x1, x2, x3, x4, R75, Rc, theta, rhov, rhol, rhop, hfg, qvolumetric, hs, tau, mdot, alphavb, dpc, epor, vp;
  double n, lref, cb, A, beta, g, jj_g;
  int i;
  // ! note: alpha=0.7, dpc=1.5m, tau=M/GA, qvolumetric=q_db_tot/V_db_tot
  n=0.5;
  lref=2.3;
  g = grav;
  cb = 0.351;
  A = 19.267;
  beta = 1.061;
  md.db.Vtot = 0.0;
  // do i = 1, md.db.ndr
  for (i=1; i<=md.db.ndr; i++){
    if(md.db.dbr[i-1].ihtmod[0] > 0.0) {continue;}
    md.db.Vtot = md.db.Vtot+md.db.dar[i-1]*md.db.dbr[i-1].zb[md.db.dbr[i-1].ndz];
  }
  rhov = md.clp.rhov;    // ! gas density
  rhol = md.clp.rhol;    // ! liquid density
  rhop = md.mp.mpp.rhom; // ! particle density
  hfg = md.clp.hfg;      // ! latent heat
  // !md.db.qrelbtm = md.ml.qrel + md.db.qrel
  qvolumetric = md.db.qrelbtm/md.db.Vtot;       // ! volumetric heat release rate (decay heat?)
  hs = md.plg.hhp-md.mj.llj*(1-sqrt(0.5));   // ! particel sedimentation height (postulation: area weighted height)
  tau = md.db.tau;              // ! melt release time duration[s]
  mdot = md.db.mdot;  // !corium mass flow rate 
  epor = md.db.epor; // ! porosity
  vp = md.db.vp;     // ! particle falling velocity (for each time step, vp = [(sum of all particle's kinetic energy)/(0.5*sum of all particle's mass)]^0.5 ) 
  // !------------------------------------------uncertainty variable---------------------------------// !
  jj_g = md.db.jj_g;
  alphavb = epor*jj_g; // !Void fraction * bubble rise velocity at bubble column
  dpc = md.db.cdpc*md.mj.di;
  if(dpc>2*pow((md.plg.csar/pi),0.5)) {
    dpc = 2*pow((md.plg.csar/pi),0.5);
  }
    // !--------------------------------------------------------------------------------------------------// !
  x1 = pow((rhol-rhov),2)/(rhop*rhov*hfg);
  x2 = qvolumetric*pow(hs,2)*tau/mdot;
  x3 = alphavb*pow(dpc,4)/(1-epor)/pow(vp,4);
  Rc = 0.414/0.674*pow((x1*x2*x3),0.3333);
  
  x1 = rhov*hfg/pow((rhol-rhov),2);
  x2 = mdot,2/qvolumetric/pow(hs,2);
  x3 = pow(vp,4)/alphavb/pow(dpc,4);
  theta = atan(4.127*x1*x2*x3);
  
  if(theta>pi/4) {
    theta = pi/4; // !Critical side slope angle
  }
  else if(theta<pi/18) {
    theta = pi/18; // !Lower boundary condition
  }
  // end if
  // !if(Rc<(md.plg.csar/3.141592),0.5 .and. Rc>0) {
  // !  md.plg.rr = Rc
  // !else
  // !  md.plg.rr = (md.plg.csar/3.141592),0.5
  // !end if
  md.db.angrep = theta*180/pi;
  md.db.tangrep = tan(theta);
  md.db.x1 = mdot;
  md.db.x2 = qvolumetric;
  md.db.x3 = hs;
  md.db.x4 = vp;
  md.db.x5 = alphavb;
  // !md.db.x6 = 
  // !md.db.x7 = 
  // !md.db.x8 = 
}




































/*###############################################################################*/ 
void mjetevol(meltdebr &md, double time, double dt)//; 
{
    double zeta, dzj, timei, timep, arj, ment, meflx, xran, x;
    double lbr1, lbr2, lbr3, rv, drdpv, drdtv, hv, dhdpv,dhdtv, rl, drdpl, drdtl, 
            hl, dhdpl, dhdtl, jacob_v, jacob_l, q_jacob, expo_jacob, 
            kepartot=0.0, kepar, parmasstot=0.0, x1, x2, x3, x4, vv;
    int i,nj, ioutp;
    int irngv, irngl;
    static int outpcnt = 0;

    // jet temperature update considering decay heat generation - juwook 220414
    md.mj.mpp.enem = md.mj.mpp.enem + qdecaymass*dt;
    menetemp(md.mppb, md.mj.mpp.enem, md.mj.temp0);

    // 216 ! check if jet flows at all // modified & added d0 <= 0 by Juwook 220414
    if (md.mj.d0*md.mj.v0 <= 1e-3){
        md.mj.zjl = md.plg.hh0;
        md.mj.vi=0, md.mj.di=0, md.mj.db = 0, md.mj.llj = 0, md.mj.llbr = 0;
        // printf("2305-mjetevol\n");
        return;
    }
    if (md.mj.d0*md.mj.d0 <= 0){
        md.mj.zjl = md.plg.hh0;
        md.mj.vi=0, md.mj.di=0, md.mj.db = 0, md.mj.llj = 0, md.mj.llbr = 0;
        // printf("2305-mjetevol\n");
        return;
    }

    // 223 ! check water level and find leading edge position << by Juwook 2022
    ttabeval(md.vjin, 0.0, vv);
    timei = -vv/grav + sqrt(pow(vv,2)+2.0*grav*(md.plg.hh0-md.plg.hhp0))/grav;
    // timei = -md.mj.v0/grav + sqrt(md.mj.v0*2 + 2*grav*md.plg.hhf)/grav; // time of water surface contact
    if (time<timei){
        md.mj.zjl = md.plg.hh0 - (md.mj.v0 + 0.5*grav*time)*time;
        md.mj.vi=0, md.mj.di=0, md.mj.db = 0, md.mj.llj = 0, md.mj.llbr = 0;
        return; 
    }

    // 231 ! now, jet is under water
    timep = time-timei; // penetration time
    md.mj.vi = md.mj.v0 + grav * timei; // velocity at water surface
    md.mj.di = md.mj.d0 * sqrt(md.mj.v0/md.mj.vi); // dia. there
    // //  mj.di = mj.d0/(1d0+2d0*grav*pg.hhf/mj.v0**2)**0.25 // dia. there (this form makes larger error)

    // 237 ! temperature assumed not change, get melt properties for jet << by Juwook 2022
    // md.mj.tempi = md.mj.temp0;
    // md.mj.mpp.tm = md.mj.tempi; // to includ 'mpp set'.. use dictionary type...
    // getmprop(md.mppb, md.mj.mpp);

    // 242 ! breakup length << Juwook 2022
    wrsteamtab(md.clp.ptot,md.clp.tspt,1,rv,drdpv,drdtv,hv,dhdpv,dhdtv,irngv);
    wrsteamtab(md.clp.ptot,md.clp.tl,0,rl,drdpl,drdtl,hl,dhdpl,dhdtl,irngl);
    jacob_v = dhdtv*(md.mj.temp0-md.clp.tspt)/md.clp.hfg;
    jacob_l = dhdtl*(md.mj.temp0-md.clp.tl)/md.clp.hfg;
    q_jacob = jacob_v/jacob_l;
    expo_jacob = 0.5*(1.0-1.0/(21354.0*q_jacob+1.0));
    lbr1 = 10.0 * sqrt(md.mj.mpp.rhom/md.clp.rhol) * md.mj.di;
    lbr2 = 2.1 * sqrt(md.mj.mpp.rhom/md.clp.rhol/grav/md.mj.di)*md.mj.vi * md.mj.di;
    lbr3 = 3.3*md.mj.di*sqrt(md.mj.mpp.rhom/md.clp.rhol)*pow((pow(md.mj.vi,2)/grav/md.mj.di),expo_jacob);

    if(md.mj.ibr == 1) // ! Taylor/Epstein type correlation with C=10. (Moriyama et al.,2005)
        {md.mj.llbr = lbr1;} // 10 * sqrt(md.mj.mpp.rhom/md.clp.rhol) * md.mj.di;}
    else if (md.mj.ibr == 2)  // ! Saito et al. (1988) correlation
        {md.mj.llbr = lbr2;} //2.1 * sqrt(md.mj.mpp.rhom/md.clp.rhol/grav/md.mj.di)*md.mj.vi * md.mj.di;}
    else if (md.mj.ibr == 3)  // ! Enhanced jet break-up model (MATE) << by Juwook 2022
        {md.mj.llbr = lbr3;}
    else if (md.mj.ibr == 4)  // ! Lbr = min (Epstein, Saito, Mate) << by Juwook 2022
        {md.mj.llbr = max(max(lbr1,lbr2),lbr3);}
    else {
        printf("mjetevol: unknown option for jet breakup\n");
        exit(EXIT_FAILURE);
    }
    md.mj.llbr = md.mj.cbr * md.mj.llbr;

    // 254 ! jet length under water, leading edge position
    md.mj.llj = min(min(md.mj.vi*timep, md.mj.llbr), md.plg.hhp);
    md.mj.zjl = md.plg.hhp - md.mj.llj;

    // 258 ! reaching bottom ? yes => deposit mass to melt lump 
    if(md.mj.llbr > md.plg.hhp && md.mj.llj >= md.plg.hhp)
    {
        md.mj.vb = md.mj.vi;
        md.mj.db = (1.0 - md.plg.hhp/md.mj.llbr) * md.mj.di;  // ! jet dia. at pool bottom
        depomlumpjet(md.mj, md.ml, md.mppb, dt);
    }
    else
    {
        md.mj.vb = 0.0; 
        md.mj.db = 0.0;
    }
            
    // 267 ! jet profile, mass entrainment, particle size << by Juwook 2022
    // ! if llj is short, use dzjet+ and calculate nj, otherwise use njet and calculate dzj
    // ! jet is discretized into nj nodes of length dzj
    nj = md.mj.njet;
    dzj = md.mj.llj / nj;
    // if(md.mj.llj < md.mj.dzjet * md.mj.njet)
    // {
    //     nj = floor( md.mj.llj / md.mj.dzjet);
    //     if (nj == 0) 
    //     {
    //         nj = 1;
    //     }
    //     dzj = md.mj.llj / nj;
    // }
    // else
    // {
    //     nj = md.mj.njet;
    //     dzj = md.mj.llj / nj;
    // }
    // end if

    // 280 ! mass entrainment flux (kg/m2s), constant over whole jet surface << by Juwook 2022
    meflx = 0.5*md.mj.mpp.rhom*md.mj.di*md.mj.vi/md.mj.llbr;
    if (md.mj.llj < md.plg.hhp) {
        meflx = 0.5*md.mj.mpp.rhom*md.mj.di*md.mj.vi/sqrt(pow(md.mj.llj,2) + 0.25*pow(md.mj.di,2));
    }
    else {
        meflx = 0.5*md.mj.mpp.rhom*md.mj.di*md.mj.vi/sqrt(pow(md.mj.llbr,2) + 0.25*pow(md.mj.di,2));
    }

    // 283 ! particle output time check & control (needed for output in trackpar)

    if(time >= md.toutp && md.dtoutp > 0.0 && outpcnt <= md.noutp){
        ioutp = 1; // ! do the output
    }
    else{
        ioutp = 0; // ! no
    }
    // end if

    // 290 ! particle generation at sections of jet and tracking of them
    // DO LOOP
        // SR: parsizegen (md, dmmsv), random_number(xran), trackpar(md,time,ioutp)
    // i = 1;
    // while i < nj:
    // << by Juwook 2022
    md.db.mdot = 0.0;
    for (i = 1; i <= nj; i++)
    {
        x1=0.0;
        x2=0.0;
        x3=0.0;
        x4=0.0;
        // ! zeta: position on a jet node, measured from water surface
        // ! order: bottom to top of the jet (close to floor first)
        zeta = md.plg.hhp - (md.mj.zjl + dzj*(i-0.5));

        // << by Juwook 2022
        arj = pi*(1.0 - zeta/md.mj.llbr)*md.mj.di * dzj; // ! side surface area (projected to vertical cylinder) (m2)
        // if(md.mj.llj < md.plg.hhp) {
        //     // arj = pi*(1d0 - zeta/md.mj.llj)*md.mj.di * dzj
        //     // arj = pi*(i-0.5)/nj*md.mj.di * dzj
        //     x1 = i*sqrt(pow((md.mj.llj*i/nj),2)+0.25*pow(md.mj.di*i/nj,2));
        //     x2 = (i-1)*sqrt(pow((md.mj.llj*(i-1)/nj),2)+0.25*pow(md.mj.di*(i-1)/nj,2));
        //     arj = 0.5*pi*md.mj.di/nj*(x1-x2);
        // }
        // else {
        //     x1 = md.mj.di*((1-md.plg.hhp/md.mj.llbr)*(1-i/nj)+i/nj);
        //     x3 = x1*sqrt(pow((md.mj.llbr+md.plg.hhp*(i/nj-1)),2)+0.25*pow(x1,2));
        //     x2 = md.mj.di*((1-md.plg.hhp/md.mj.llbr)*(1-(i-1)/nj)+(i-1)/nj);
        //     x4 = x2*sqrt(pow((md.mj.llbr+md.plg.hhp*((i-1)/nj-1)),2)+0.25*pow(x2,2));

        //     arj = 0.5*pi*(x3-x4);
        //     // arj = pi*(1d0 - zeta/md.mj.llbr)*md.mj.di * dzj
        // }
        // // endif

        ment = meflx * arj * dt; // ! mass entrainment (kg) from position zeta (length dzj)
        md.db.mdot = md.db.mdot + meflx * arj; // << Juwook 2022
        // 298 ! melt particle size and other properties
        parsizegen(md, md.mp.dmmb);
        md.mp.mpp = md.mj.mpp;
        md.mp.mass = pi/6.0*pow(md.mp.dia,3)*md.mp.mpp.rhom;
        md.mp.ars = pi*pow(md.mp.dia,2);
        md.mp.num = ment/md.mp.mass;
        md.mp.z = md.plg.hhp - zeta; 
        md.mp.r = 0.0;
        // 305 ! initial velocity of particles
        md.mp.vz = md.mp.cvz * md.mj.vi;

        if(md.mp.iranvp == 1){
            xran = randomX(); // call random_number(xran)
            md.mp.vz = md.mp.vz * xran;
        } // end if
        md.mp.vr = md.mp.cvr * md.mp.vz;
        if(md.mp.iranvp == 1){
            xran = randomX(); //call random_number(xran)
            md.mp.vr = md.mp.vr * xran;
        } // end if
        md.mp.vz = - md.mp.vz; // ! z ordinate is directing upward !!

        md.masspartot = md.masspartot + md.mp.mass*md.mp.num;

        // ! track the generated particle to evaluate contribution as transient heat source
        // ! and debris bed formation
        // ! (need vapor properties for film boiling model, i.e. film props. rhovf,muvf,cpvf,lamvf)
        
        md.mp.id = md.mp.id + 1; // ! give a number to a new particle

        trackpar(md,time,ioutp,kepar); // << by Juwook 2022
        // call trackpar(md,time,ioutp)  // ! particle data output done in trackpar
                                        // ! this may cause output of aborted data
                                        // ! in case of failure in solution of ho
        parmasstot = parmasstot + md.mp.mass*md.mp.num;
        kepartot = kepartot + md.mp.num*0.5*md.mp.mass*kepar;
    }
    // end do

    // << by Juwook 2022
    md.db.vp = pow((kepartot/0.5/parmasstot),0.5);
    if(md.db.ishape==2 && md.db.qrelbtm>0 && md.db.jj_g>0) {
        conical(md,dt);
    }
    // end if

    // 330 ! control particle data output time after output is done
    if(ioutp > 0){
        md.toutp = md.toutp + md.dtoutp;
        outpcnt = outpcnt + 1;
    }
    // end if
    // 336 ! tidy debris bed (set vertical node boundary positions for calculated node sizes)
        // SR fixdbednode (md%db)
    fixdbednode(md.db);

    // 339 ! total jet inlet: mass and energy with base temperature 300K (enem0base)
    x = pi/4.0*pow(md.mj.d0,2) * md.mj.v0 * md.mj.mpp.rhom * dt;
    md.mj.massintot = md.mj.massintot + x;
    md.mj.eneintot = md.mj.eneintot + x*max(md.mj.mpp.enem - enem0base, 0.0);
    return; 
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void initmlump(meltlmp &ml)
{
    ml.mass = 0.0;
    ml.thk = 0.0;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void initmhex (melthex &mh)
{
    int i;
    mh.thex[0]=0.0;
    // mh.qhex[0] = 0.0;
    for (i=1; i<=nthexmax; i++)
    {
        mh.thex[i] = mh.thex[i-1]+mh.dthex;
        mh.qhex[i-1] = 0.0;
    }
    mh.nthex = 0; 

}
/*###############################################################################*/ 

/*###############################################################################*/ 
void decayheat(double time)
//   ! only the calculation time is given as an argument
//   ! following parameters and outputs are shared as global vars.
//   ! time after shutdown at time=0 (s): t0sdown,
//   ! operation power (W) : powopr, total core debris mass (kg): coremass,
//   ! decay heat rate by fraction to operation power : qdcfrac, 
//   ! total decay heat power (W): qdecay, 
//   ! decay heat fraction outside the coremass : fqdcex
//   ! decay heat power per mass of core debris (W/kg): qdecaymass
{
    qdcfrac = 0.1250 * pow((time + t0sdown),(-0.2752));
    qdecay = qdcfrac * powopr;
    if (coremass > 0) {
        qdecaymass = qdecay/coremass*(1-fqdcex); // ! power per mass of core debris
    }
    else{                                        // ! qdecay-qdecaymass*(debris+lump mass)
        qdecaymass = 0;                          // ! should be counted inside containment
    }
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void initdecayheat(double t0sd,double pow,double cmass,double fex)
// ! initialization (paremeter setting) for decay heat model
{
    // ! time after shutdown at time=0 (s), operation power (W), core mass (kg)
    t0sdown = t0sd;
    powopr = pow;
    coremass = cmass;
    fqdcex = fex;
}
/*###############################################################################*/ 


















/*###############################################################################*/ 
void dbstackcollapse(dbstack &dbs, int i)
{
    int nz, k;


    nz = dbs.ndz; 
    if(i<1 || i>nz) {
        // call errorhdr(mess="dbstackcollapse: invalid index",act="stop")
        printf("dbstackcollapse: invalid index\n");
        exit(EXIT_FAILURE);
    }
    if(i < nz) {
        for (k=i; k<=nz-1; k++) // do k = i, nz-1
        {
            dbs.dz[k-1] = dbs.dz[k];
            dbs.por[k-1] = dbs.por[k];
            dbs.mass[k-1] = dbs.mass[k];
            dbs.num[k-1] = dbs.num[k];
            dbs.ars[k-1] = dbs.ars[k];
            dbs.dia[k-1] = dbs.dia[k];
            dbs.tsf[k-1] = dbs.tsf[k];
            dbs.qrel[k-1] = dbs.qrel[k];
            dbs.ihtmod[k-1] = dbs.ihtmod[k];
            dbs.mpp[k-1] = dbs.mpp[k];
        }
    }
    dbs.zb[nz] = dbs.zb[nz-1];
    dbs.zc[nz-1] = dbs.zb[nz-1];
    dbs.dz[nz-1] = 0.0;
    dbs.mass[nz-1] = 0.0;
    dbs.num[nz-1] = 0.0;
    dbs.ars[nz-1] = 0.0;
    dbs.tsf[nz-1] = 0.0;
    dbs.qrel[nz-1] = 0.0;
    dbs.mpp[nz-1].enem = 0.0;
    dbs.ihtmod[nz-1] = 0;
    dbs.ndz = max(nz - 1, 1);
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void dbrepose(debrbed &db, mpropbase mppb)
{
    double xtan, ar0,ar1,s0,s1,dzm,dz,x0,x1;
    int i, nr, iavl, n0, n1;

    nr = db.ndr;
    do
    {
        iavl = 0; // ! flag for occurrence of avalanche
        for (i=1; i<=nr-1; i++) // do i = 1, nr-1
        {
            // ! tangent of the slop i..i+1
            xtan = (db.dbr[i-1].zb[db.dbr[i-1].ndz] - db.dbr[i].zb[db.dbr[i].ndz]) 
                /(db.rc[i]-db.rc[i-1]);
            if(xtan > db.tangrep) { // ! avalanche occurs
                iavl = 1;
                // ! move node i,ndz to i+1
                n0 = db.dbr[i-1].ndz;
                n1 = db.dbr[i].ndz;
                ar0=db.dar[i-1]; // ! pi*(db.rb[i-1]**2 - db.rb(i-1)**2) // ! base area of origin
                ar1=db.dar[i]; // ! pi*(db.rb[i]**2 - db.rb[i-1]**2) // ! base area of destination
                dzm = ar0/ar1*db.dbr[i-1].dz[n0-1]; // ! height of the top node after moved
                dz = db.dbr[i].dz[n1-1];       // ! height of the top node at i+1
                if(dzm + dz > db.dzdb) { // ! new dz[i] is big, a new node is generated
                    db.dbr[i].ndz = n1 + 1;
                    db.dbr[i].dz[n1]     =  dzm;
                    db.dbr[i].por[n1]    =  db.dbr[i-1].por[n0-1];
                    db.dbr[i].mass[n1]   =  db.dbr[i-1].mass[n0-1];
                    db.dbr[i].num[n1]    =  db.dbr[i-1].num[n0-1];
                    db.dbr[i].ars[n1]    =  db.dbr[i-1].ars[n0-1];
                    db.dbr[i].dia[n1]    =  db.dbr[i-1].dia[n0-1];
                    db.dbr[i].tsf[n1]    =  db.dbr[i-1].tsf[n0-1];
                    db.dbr[i].qrel[n1]   =  db.dbr[i-1].qrel[n0-1];
                    db.dbr[i].ihtmod[n1] =  db.dbr[i-1].ihtmod[n0-1];
                    db.dbr[i].mpp[n1]    =  db.dbr[i-1].mpp[n0-1];
                }
                else                               // ! not big, merge them as one node
                {
                    // ! pore volume
                    x0 = db.dbr[i].por[n1-1]*db.dbr[i].dz[n1-1]*ar1 + db.dbr[i-1].por[n0-1]*db.dbr[i-1].dz[n0-1]*ar0;
                    x1 = db.dbr[i].dz[n1-1]*ar1 + db.dbr[i-1].dz[n0-1]*ar0;
                    db.dbr[i].por[n1-1] = x0/x1;
                    // ! surface area weighted av. for tsf
                    s0 = db.dbr[i-1].ars[n0-1];
                    s1 = db.dbr[i].ars[n1-1];
                    db.dbr[i].tsf[n1-1] = (s1*db.dbr[i].tsf[n1-1] + s0*db.dbr[i-1].tsf[n0-1])/(s1+s0);
                    // ! qrel
                    db.dbr[i].qrel[n1-1] = db.dbr[i].qrel[n1-1] + db.dbr[i-1].qrel[n0-1];
                    // ! ihtmod: retain if same, otherwise clear
                    if(db.dbr[i].ihtmod[n1-1] != db.dbr[i-1].ihtmod[n0-1]) {db.dbr[i].ihtmod[n1-1] = 0;}

                    // ! internal energy and temperature
                    x0 = (db.dbr[i].mpp[n1-1].enem*db.dbr[i].mass[n1-1] + db.dbr[i-1].mpp[n0-1].enem* 
                        db.dbr[i-1].mass[n0-1])/(db.dbr[i].mass[n1-1] + db.dbr[i-1].mass[n0-1]);
                    menetemp(mppb,x0,x1);
                    db.dbr[i].mpp[n1-1].tm = x1;
                    // ! debris phys. props. at updated temperature
                    getmprop(mppb,db.dbr[i].mpp[n1-1]);

                    // ! volume at new temperature
                    x0 = db.dbr[i-1].mass[n0-1]/db.dbr[i].mpp[n1-1].rhom;
                    x1 = db.dbr[i].mass[n1-1]/db.dbr[i].mpp[n1-1].rhom;
                    // ! mass update
                    db.dbr[i].mass[n1-1] = db.dbr[i].mass[n1-1] + db.dbr[i-1].mass[n0-1];
                    // ! surface area update
                    s0 = pow((36.0*pi*pow(x0,2)*db.dbr[i-1].num[n0-1]),(1.0/3.0));
                    s1 = pow((36.0*pi*pow(x1,2)*db.dbr[i].num[n1-1]),(1.0/3.0));
                    db.dbr[i].ars[n1-1] = s0+s1;
                    // ! number update
                    db.dbr[i].num[n1-1] = pow((s0+s1),3)/pow((x0+x1),2)/36.0/pi;
                    // ! particle diameter update
                    db.dbr[i].dia[n1-1] = 6.0*(x0+x1)/(s0+s1);
                }
                dbstackcollapse(db.dbr[i-1],db.dbr[i-1].ndz);
                fixdbednode(db);
            }
        }
    }
    while (iavl != 0);



    // << by Juwook 2022
    // reverse avalanche
    do{
        iavl = 0; // ! flag for occurrence of avalanche
        for (i=1; i<=nr-1; i++) // do i = 1, nr-1
        {
            // ! tangent of the slop i..i+1
            xtan = (db.dbr[i].zb[db.dbr[i].ndz] - db.dbr[i-1].zb[db.dbr[i-1].ndz]) 
                /(db.rc[i]-db.rc[i-1]);
            if(xtan > db.tangrep) { // ! avalanche occurs
                iavl = 1;
                // ! move node i,ndz to i+1
                n0 = db.dbr[i-1].ndz;
                n1 = db.dbr[i].ndz;
                ar0=db.dar[i-1]; // ! pi*(db.rb[i-1]**2 - db.rb(i-1)**2) // ! base area of origin
                ar1=db.dar[i]; // ! pi*(db.rb[i]**2 - db.rb[i-1]**2) // ! base area of destination
                dzm = ar1/ar0*db.dbr[i].dz[n1-1]; // ! height of the top node after moved
                dz = db.dbr[i-1].dz[n0-1];       // ! height of the top node at i+1
                if(dzm + dz > db.dzdb) { // ! new dz[i] is big, a new node is generated
                    db.dbr[i-1].ndz = n0 + 1;
                    db.dbr[i-1].dz[n0]     =  dzm;
                    db.dbr[i-1].por[n0]    =  db.dbr[i].por[n1-1];
                    db.dbr[i-1].mass[n0]   =  db.dbr[i].mass[n1-1];
                    db.dbr[i-1].num[n0]    =  db.dbr[i].num[n1-1];
                    db.dbr[i-1].ars[n0]    =  db.dbr[i].ars[n1-1];
                    db.dbr[i-1].dia[n0]    =  db.dbr[i].dia[n1-1];
                    db.dbr[i-1].tsf[n0]    =  db.dbr[i].tsf[n1-1];
                    db.dbr[i-1].qrel[n0]   =  db.dbr[i].qrel[n1-1];
                    db.dbr[i-1].ihtmod[n0] =  db.dbr[i].ihtmod[n1-1];
                    db.dbr[i-1].mpp[n0]    =  db.dbr[i].mpp[n1-1];
                }
                else                               // ! not big, merge them as one node
                {
                    // ! pore volume
                    x0 = db.dbr[i-1].por[n0-1]*db.dbr[i-1].dz[n0-1]*ar0 + db.dbr[i].por[n1-1]*db.dbr[i].dz[n1-1]*ar1;
                    x1 = db.dbr[i-1].dz[n0-1]*ar0 + db.dbr[i].dz[n1-1]*ar1;
                    db.dbr[i-1].por[n0-1] = x0/x1;
                    // ! surface area weighted av. for tsf
                    s0 = db.dbr[i-1].ars[n0-1];
                    s1 = db.dbr[i].ars[n1-1];
                    db.dbr[i-1].tsf[n0-1] = (s0*db.dbr[i-1].tsf[n0-1] + s1*db.dbr[i].tsf[n1-1])/(s1+s0);
                    // ! qrel
                    db.dbr[i-1].qrel[n0-1] = db.dbr[i-1].qrel[n0-1] + db.dbr[i].qrel[n1-1];
                    // ! ihtmod: retain if same, otherwise clear
                    if(db.dbr[i-1].ihtmod[n0-1] != db.dbr[i].ihtmod[n1-1]) {db.dbr[i-1].ihtmod[n0-1] = 0;}

                    // ! internal energy and temperature
                    x0 = (db.dbr[i-1].mpp[n0-1].enem*db.dbr[i-1].mass[n0-1] + db.dbr[i].mpp[n1-1].enem* 
                        db.dbr[i].mass[n1-1])/(db.dbr[i].mass[n1-1] + db.dbr[i-1].mass[n0-1]);
                    menetemp(mppb,x0,x1);
                    db.dbr[i-1].mpp[n0-1].tm = x1;
                    // ! debris phys. props. at updated temperature
                    getmprop(mppb,db.dbr[i-1].mpp[n0-1]);

                    // ! volume at new temperature
                    x0 = db.dbr[i-1].mass[n0-1]/db.dbr[i-1].mpp[n0-1].rhom;
                    x1 = db.dbr[i].mass[n1-1]/db.dbr[i-1].mpp[n0-1].rhom;
                    // ! mass update
                    db.dbr[i-1].mass[n0-1] = db.dbr[i-1].mass[n0-1] + db.dbr[i].mass[n1-1];
                    // ! surface area update
                    s0 = pow((36.0*pi*pow(x0,2)*db.dbr[i-1].num[n0-1]),(1.0/3.0));
                    s1 = pow((36.0*pi*pow(x1,2)*db.dbr[i].num[n1-1]),(1.0/3.0));
                    db.dbr[i-1].ars[n0-1] = s0+s1;
                    // ! number update
                    db.dbr[i-1].num[n0-1] = pow((s0+s1),3)/pow((x0+x1),2)/36.0/pi;
                    // ! particle diameter update
                    db.dbr[i-1].dia[n0-1] = 6.0*(x0+x1)/(s0+s1);
                }
                dbstackcollapse(db.dbr[i],db.dbr[i].ndz);
                fixdbednode(db);
            }
        }
    }
    while (iavl != 0);
        // if (iavl == 0) exit




}
/*###############################################################################*/ 

/*###############################################################################*/ 
void dbrepose_in(debrbed &db, mpropbase mppb)
{
    // ! debris bed avalanche into the in-side direction (optional application)
    double xtan, ar0,ar1,s0,s1,dzm,dz,x0,x1;
    int i, nr, iavl, n0, n1;


    nr = db.ndr;
    do
    {
        iavl = 0; // ! flag for occurrence of avalanche
        for (i=nr; i==2; i--) // do i = nr, 2, -1
        {
            // ! tangent of the slop i..i+1
            xtan = (db.dbr[i-1].zb[db.dbr[i-1].ndz] - db.dbr[i-2].zb[db.dbr[i-2].ndz]) 
                / (db.rc[i-1]-db.rc[i-2]);
            // ! move node i,ndz to i-1
            // ! predict height of receiving ring after moving and if the slope becomes 
            // ! opposite, do not move (in-side avalanche can make too high destination)
            n0 = db.dbr[i-1].ndz;
            n1 = db.dbr[i-2].ndz;
            ar0=db.dar[i-1]; // ! pi*(db.rb[i-1]**2 - db.rb[i]**2) // ! base area of origin
            ar1=db.dar[i-2]; // !pi*(db.rb[i]**2 - db.rb(i-2)**2) // ! base area of destination
            dzm = ar0/ar1*db.dbr[i-1].dz[n0-1]; // ! height of the top node after moved
            dz = db.dbr[i-2].dz[n1-1];        // ! height of the top node at i-1
            x0 = db.dbr[i-1].zb[n0-1] - (db.dbr[i-2].zb[n1] + dzm); // ! gap after move
            if(xtan > db.tangrepin && x0 > 0.0) { // ! avalanche occurs
                iavl = 1;
                if(dzm + dz > db.dzdb) { // ! new dz[i] is big, a new node is generated
                    db.dbr[i-2].ndz = n1 + 1;
                    db.dbr[i-2].dz[n1]     =  dzm;
                    db.dbr[i-2].por[n1]    =  db.dbr[i-1].por[n0-1];
                    db.dbr[i-2].mass[n1]   =  db.dbr[i-1].mass[n0-1];
                    db.dbr[i-2].num[n1]    =  db.dbr[i-1].num[n0-1];
                    db.dbr[i-2].ars[n1]    =  db.dbr[i-1].ars[n0-1];
                    db.dbr[i-2].dia[n1]    =  db.dbr[i-1].dia[n0-1];
                    db.dbr[i-2].tsf[n1]    =  db.dbr[i-1].tsf[n0-1];
                    db.dbr[i-2].qrel[n1]   =  db.dbr[i-1].qrel[n0-1];
                    db.dbr[i-2].ihtmod[n1] =  db.dbr[i-1].ihtmod[n0-1];
                    db.dbr[i-2].mpp[n1]    =  db.dbr[i-1].mpp[n0-1];
                }
                else {                       // ! not big, merge them as one node
                    // ! pore volume
                    x0 = db.dbr[i-2].por[n1-1]*db.dbr[i-2].dz[n1-1]*ar1 + db.dbr[i-1].por[n0-1]*db.dbr[i-1].dz[n0-1]*ar0;
                    x1 = db.dbr[i-2].dz[n1-1]*ar1 + db.dbr[i-1].dz[n0-1]*ar0;
                    db.dbr[i-2].por[n1-1] = x0/x1;
                    // ! surface area weighted av. for tsf
                    s0 = db.dbr[i-1].ars[n0-1];
                    s1 = db.dbr[i-2].ars[n1-1];
                    db.dbr[i-2].tsf[n1-1] = (s1*db.dbr[i-2].tsf[n1-1] + s0*db.dbr[i-1].tsf[n0-1])/(s1+s0);
                    // ! qrel
                    db.dbr[i-2].qrel[n1-1] = db.dbr[i-2].qrel[n1-1] + db.dbr[i-1].qrel[n0-1];
                    // ! ihtmod: retain if same, otherwise clear
                    if(db.dbr[i-2].ihtmod[n1-1] != db.dbr[i-1].ihtmod[n0-1]) {db.dbr[i-2].ihtmod[n1-1] = 0;}

                    // ! internal energy and temperature
                    x0 = (db.dbr[i-2].mpp[n1-1].enem*db.dbr[i-2].mass[n1-1] + db.dbr[i-1].mpp[n0-1].enem* 
                        db.dbr[i-1].mass[n0-1])/(db.dbr[i-2].mass[n1-1] + db.dbr[i-1].mass[n0-1]);
                    menetemp(mppb,x0,x1);
                    db.dbr[i-2].mpp[n1-1].tm = x1;
                    // ! debris phys. props. at updated temperature
                    getmprop(mppb,db.dbr[i-2].mpp[n1-1]);

                    // ! volume at new temperature
                    x0 = db.dbr[i-1].mass[n0-1]/db.dbr[i-2].mpp[n1-1].rhom;
                    x1 = db.dbr[i-2].mass[n1-1]/db.dbr[i-2].mpp[n1-1].rhom;
                    // ! mass update
                    db.dbr[i-2].mass[n1-1] = db.dbr[i-2].mass[n1-1] + db.dbr[i-1].mass[n0-1];
                    // ! surface area update
                    s0 = pow((36.0*pi*pow(x0,2)*db.dbr[i-1].num[n0-1]),(1.0/3.0));
                    s1 = pow((36.0*pi*pow(x1,2)*db.dbr[i-2].num[n1-1]),(1.0/3.0));
                    db.dbr[i-2].ars[n1-1] = s0+s1;
                    // ! number update
                    db.dbr[i-2].num[n1-1] = pow((s0+s1),3)/pow((x0+x1),2)/36.0/pi;
                    // ! particle diameter update
                    db.dbr[i-2].dia[n1-1] = 6.0*(x0+x1)/(s0+s1);
                }
                dbstackcollapse(db.dbr[i-1],db.dbr[i-1].ndz);
                fixdbednode(db);
            }
        }
    }
    while (iavl !=0);
    // if (iavl == 0) exit
    // end do

}
/*###############################################################################*/ 



/*###############################################################################*/ 
void dbrlmpevol(meltdebr &md, double time, double dt)
{
    int i, j, ndr, ndz, imod, imodel;
    double qhtr, qdhf, qdbr, pr,tsat,tsat1,tvf, 
        rhovf,cpvf,lamvf,muvf, qflx, tav,tavn,etav,tliq,etliq,tcr,etcr,dedt,
        mrm,trm,erm,qrel,tsf, x0,x1,x2,x3,x4, dbrad, qbtm, qbtmtot, alpini = 0.0,
        mrmtot;
    double qdbrtot, Atot, Vtot, jj_g, Adhf;
    // ! dhf model control vars: db%idhf, db%alpinidhf

    if(md.db.idhf == 1){imodel = 3;} // ! Rahman model
    else if(md.db.idhf ==2){imodel = 1;} // ! Reed model
    else if(md.db.idhf ==3){imodel = 4;} // ! MLee model (Mooneon Lee 2017)
    else if(md.db.idhf != 0){
        // call errorhdr(mess='dbrlmpevol: unknown DHF model',act='stop')
        printf("dbrlmpevol: unknown DHF model\n");
        exit(EXIT_FAILURE);
    }
    if(md.db.alpinidhf > 0.0) {alpini = md.db.alpinidhf;}
    
    // ! melt lump temperature evolusion: phase 1
    // ! this phase does not include heat transfer from the lump top surface
    // ! (needed before addition from debris bed remelting where 
    // !  decay heat is included already)
    x0 = qdecaymass*dt;
    md.ml.mpp.enem = md.ml.mpp.enem + x0;
    menetemp(md.mppb,md.ml.mpp.enem,md.ml.mpp.tm);
    getmprop(md.mppb,md.ml.mpp);
    md.ml.qdctot = md.ml.qdctot + x0*md.ml.mass;

    // ! save lump decayheat (W)
    md.ml.qdc = qdecaymass*md.ml.mass;

    // ! lump thickness (area is assumed the same as existing debris bed)
    dbrad = 0.0; // ! radius of debris bed, assume lump has the same size
    md.ml.thk = 0.0; // ! clear thickness
    for (i=1; i<=md.db.ndr; i++)//    do i = 1,md.db.ndr
    {
        if(md.db.dbr[i-1].mass[0] > 1e-3) {dbrad = md.db.rb[i];}
    }
    if(dbrad > 1e-3) {md.ml.thk = md.ml.mass/md.ml.mpp.rhom/(pi*pow(dbrad,2));}
    md.ml.thk = max(md.ml.thk,1e-2); // ! give 1cm as a minimum heat transfer thickness

    // ! 1D DHF evaluatoion vs particle heat transfer for every ring
    // !
    // ! debris bed data structure/contents
    // !  i = 1 .. md.db.ndr : index of rings
    // !  md.db.dbr(1:md.db.ndr) : a vertical column of debris elements
    // !  j = 1 .. md.db.dbr(i).ndz : index of vertical elements in the ring node i
    // !  md.db.dbr(i).por(j), mass(j), num(j), ars(j), dia(j) : porosity, mass, 
    // !   number of particles, surface area of debris in the node i,j
    // !  md.db.dbr(i).mpp(j) : debris physical properties in the node i,j
    // !   ..mpp(j).tm, tmel, tsol, tliq, rhom, enem, cpm, dedt, lamm, epsm
    // !    T, solidus T, liquidus T, density, internal E, spec. heat, 
    // !    spec. heat considering latent heat, heat conductivity, emmissivity
    // !    (viscosity, mum, and surface tension, sigm, are not needed in solid debris model)
    // !
    // ! for every ring node of debris bed
    // !   calc DHF based max. removable power : qdhf (W)
    // !   calc heat transfer based power (particle HT * attenuation factor) : qhtr (W)
    // !   take minimum of the above two for actual power : qdbr (W)
    // !   apply qdbr/total_mass for energy conservation of all vertical nodes
    // !   add heat release (J) to heat transfer buffer by depomhex(..)
    // !   re-arrange the mass/energy of the node to lump if remelting occurs
    // !   update debris bed stack node properties for the new temperature
    // !   fix vertical boundary positions and check for avalanche

    ndr = md.db.ndr;
    qbtmtot = 0.0; // ! sum of melt lump heat release rate
    mrmtot = 0.0; // ! sum of debris mass remelting
    // << juwook 2022
    qdbrtot = 0.0;
    Atot = 0.0;
    Adhf = 0.0;

    md.db.qrel = 0.0; 
    md.db.qdc = 0.0; // ! debris heat release and decay heat (W)
    // << juwook 2022
    md.db.jj_g = 0.0; 
    md.db.qrelbtm = 0.0;

    for(i=1; i<=ndr; i++) // do i = 1, ndr
    {
        if(md.db.dbr[i-1].mass[0] <= 0.0) {continue;}
        ndz = md.db.dbr[i-1].ndz;
        // ! calc heat transfer based power (particle HT * attenuation factor) : qhtr (W)
        // ! velocity is neglected, water temperature is assumed to be tsat for ptot,
        // ! radiation is neglected
        qhtr = 0.0;
        for (j=1; j<=ndz; j++) // do j = 1, ndz
        {
            // ! use tav; surface temperature drop not considered (usage of single particle 
            // ! heat transfer model itself is very crude, and already a strong empirical 
            // ! factor is used)
            tav = md.db.dbr[i-1].mpp[j-1].tm;
            tsf = md.db.dbr[i-1].tsf[j-1];
            pr = md.clp.ptot;
            tsat = md.clp.tspt;
            // !check tsf and give limits
            tsf = max(min(tsf,tav),tsat);
            // ! vapor film properties need evaluation here (at ptot,0.5*(tsat+tav))
            tvf = max( 0.5*(tsf+tsat), tsat );
            wrsteamtab(pr,tvf,1, rhovf,x0,x1,x3,x4,cpvf,imod);
            lamvf = tconwaterf(rhovf,tvf);
            muvf = viscwaterf(rhovf,tvf);
            // ! heat transfer
            htranspar(tsf,md.db.dbr[i-1].dia[j-1],md.db.dbr[i-1].mpp[j-1].epsm, 0.0,0.0, 
                    pr,0.0, tsat,tsat,tsat,md.clp.hfg,
                    md.clp.rhov,md.clp.rhol,md.clp.cpv,md.clp.cpl,
                    md.clp.muv,md.clp.mul,md.clp.lamv,md.clp.laml,md.clp.sig, 
                    rhovf,cpvf,lamvf,muvf, 0, imod, qflx, x0, x1, x2);
            qflx = qflx * md.db.chtc; // ! modify heat flux by given factor (should be 
                                    // ! a strong attenuation considering crowded situation)
            qhtr = qhtr + qflx * md.db.dbr[i-1].ars[j-1]; // ! sum of heat transfer based power (W)
            md.db.dbr[i-1].ihtmod[j-1] = imod; // ! save heat transfer mode
        }

        // ! calculate total and average vars for ring i for DHF evaluation etc.
        // ! average porosity and particle dia (sauter) in the ring i 
        x0 = 0.0; 
        x1 = 0.0; 
        x2 = 0.0; 
        x3 = 0.0;
        // ! ring totals
        for (j=1; j<=ndz; j++) // do j = 1, ndz
        {
            x0 = x0 + md.db.dbr[i-1].por[j-1] * md.db.dbr[i-1].dz[j-1]; // ! pore volume (normalized)
            x1 = x1 + md.db.dbr[i-1].mass[j-1]/md.db.dbr[i-1].mpp[j-1].rhom; // ! debris volume
            x2 = x2 + md.db.dbr[i-1].ars[j-1];  // ! debris surface area
            x3 = x3 + md.db.dbr[i-1].mass[j-1]; // ! debris mass
        }
        x0 = x0/md.db.dbr[i-1].zb[ndz]; // ! av. porosity of ring i
        x1 = 6.0*x1/x2; // ! sauter mean dia of ring i

        // ! consider q" at debris bed bottom (melt lump top)
        // ! q"dbr+q"btm < DHF if DHF is considered 
        x4 = md.db.dar[i-1]; // !pi*(md.db.rb[i-1]**2-md.db.rb(i-1)**2) // ! base area of the ring
        if(md.ml.mass > 1e-3 && md.db.rb[i] <= dbrad) {
            qbtm = md.ml.fhtlmp*2.0*md.ml.mpp.lamm/md.ml.thk 
                *max(md.ml.mpp.tm-md.clp.tspt, 0.0)*x4; // ! lump surface power (W) for the ring
        }
        else {qbtm = 0.0;}

        // ! calc DHF based max. removable power : qdhf (W)
        qdhf = 0.0;
        if(md.db.idhf > 0){
            // ! DHF model called with default options: Tsat at given P, Rahman model
            // ! and water inflow rate=0
            setprop_dhf(pr,tsat1,x0,x1,md.db.vlbtm);
            get_alp_dhf(x0,x1,imod,imodel,alpini); // ! get alp, dhf, model (1:Reed, 3:Rahman, 4:MLee)
            qdhf = x1*x4; // ! DHF based power (W)

            // ! possibly consider 2D enhance effect by factor, fdhf
            qdhf = md.db.fdhf * qdhf;

            qbtm = min(qbtm, qdhf*md.ml.fqbdhf); // ! => heat release from lump for this dbr ring
                                                // !   limited < DHF
            qdhf = qdhf - qbtm;                  // ! the actual DHF for debris
        }

        // ! take smaller of the qhtr and qdhf for actual power : qdbr (W)
        if(md.db.idhf > 0 && qdhf < qhtr) {
            qdbr = qdhf;
            for (int k=0; k<ndz;k++)
            {
                md.db.dbr[i-1].ihtmod[k] = -imod;
            }
            // md.db.dbr[i-1].ihtmod(1:ndz) = -imod;
            qdbrtot = qdbrtot + qdbr + qbtm;
            md.db.qrelbtm = md.db.qrelbtm + qdhf + qbtm;
            if(md.db.rb[i-1] <= dbrad) {
                md.db.jj_g = md.db.jj_g+jj_g*md.db.dar[i-1];
                Adhf = Adhf + md.db.dar[i-1];
            } // end if
        }
        else {
            qdbr = qhtr;
            qdbrtot = qdbrtot + qdbr; // << Juwook 2022
        }

        // ! save debris heat release power (W)
        md.db.qrel = md.db.qrel + qdbr;
        Atot = Atot +x4 ; // << juwook 2022
        qflx = qdbr/x2; // ! average surface heat flux for ring i (W/m2)
        qdbr = qdbr/x3; // ! heat release power per mass (W/kg) for ring i

        // ! save decay heat generation in ring i
        md.db.qdctot = md.db.qdctot + x3*qdecaymass*dt; // ! total (J)
        md.db.qdc = md.db.qdc + x3*qdecaymass;  // ! power (W)

        // ! sum of heat release rate from melt lump
        qbtmtot = qbtmtot + qbtm;

        // ! energy conservation eq. solution for every vertical stack node
        // ! temperature update uses implicit scheme with effective heat transfer rate
        // ! (-q + qdc)*dt = demdt*(tavn - tav), q=Reff*(Tavn - Tsat)
        // ! => Tavn = ((Reff*Tsat + qdc)*dt + demdt*Tav)/(Reff*dt + demdt)
        for (j=1; j<=ndz; j++) // do j = 1, ndz
        {
            tav = md.db.dbr[i-1].mpp[j-1].tm;
            etav = md.db.dbr[i-1].mpp[j-1].enem;
            tliq = md.db.dbr[i-1].mpp[j-1].tliq;
            tsf = md.db.dbr[i-1].tsf[j-1];
            // ! T criterion for remelting
            if(md.ml.itcr == 1) { // ! tav compared with tmel+..
                tcr = md.db.dbr[i-1].mpp[j-1].tmel + 0.5* 
                    md.ml.ftcr*(md.db.dbr[i-1].mpp[j-1].tliq-md.db.dbr[i-1].mpp[j-1].tsol);
            }
            else if(md.ml.itcr == 2) { // ! tsf compared with tsol+..
                tcr = md.db.dbr[i-1].mpp[j-1].tsol + 
                    md.ml.ftcr*(md.db.dbr[i-1].mpp[j-1].tliq-md.db.dbr[i-1].mpp[j-1].tsol);
            }
        // !      tsat = md.clp.tspt
            dedt = md.db.dbr[i-1].mpp[j-1].dedt;
            x0 = fabs(qdbr/max(tav-tsat,1.0)); // ! per mass effective heat transfer rate
            x1 = ((x0*tsat+qdecaymass)*dt + dedt*tav)/(x0*dt+dedt); // ! updated debris temp.
            tavn = x1;
            mrm = 0.0;
            trm = 0.0;
            erm = 0.0;

            // ! released heat (J) in node j, exact formula
            mtempene(md.mppb,tavn,x2);
            qrel = (etav-x2+qdecaymass*dt)*md.db.dbr[i-1].mass[j-1];

            // ! if updated temperature exceeds remelting criterion, eval. remelting mass
            // ! by alternative energy equation
            // ! m*((-q + qdc)*dt + Em(Tav)) = (m - mrm)*Em(Tcr) + mrm*Em(Tliq)
            // ! => mrm = (m*((-q + qdc)*dt + Em(Tav) - Em(Tcr)))/(Em(Tliq) - Em(Tcr))
            // ! if mrm>m yet another solution by forcing mrm=m
            // ! (-q + qdc)*dt + Em(Tav) =  Em(Tavn) : Tavn=the new temp. as liquid
            if(md.ml.itcr == 1) {x1 = tavn;}
            else if(md.ml.itcr == 2) {
                x1 = tsf + (tavn-tav); // ! rough esitimate of the new surface temp. 
            }
            if (md.ml.itcr > 0 && x1 > tcr) {
                mtempene(md.mppb,tcr,etcr);
                mtempene(md.mppb,tliq,etliq);
                // ! remelting mass with T=Tliq, while remaining solid is at T=Tcr
                x2 = md.db.dbr[i-1].mass[j-1] * ((-qdbr+qdecaymass)*dt+etav-etcr)/(etliq - etcr);
                tavn = tcr;
                mrm = x2;
                trm = tliq;
                erm = etliq;
                qrel = qdbr*dt*md.db.dbr[i-1].mass[j-1]; // ! released heat (J)
                // ! what if all mass is remelted? => new temperature in liquid zone
                if(x2 > md.db.dbr[i-1].mass[j-1]) {
                    x2 = md.db.dbr[i-1].mass[j-1]; // ! mass=all
                    // ! evaluate exact energy/temperature of liquid produced by remelting
                    x3 = (-qdbr+qdecaymass)*dt+etav; // ! x3: new internal energy (J/kg)
                    menetemp(md.mppb,x3,x4); // ! x4: new temperature (as liquid)
                    tavn = 0.0;
                    mrm = x2;
                    trm = x4;
                    erm = x3;
                    // ! qrel is the same as above
                }
            }

            // ! re-arrange debris bed and melt lump mass/energy in case of remelting
            if(mrm > 0.0) {
                x0 = md.db.dbr[i-1].mass[j-1]; // ! number to be proportional to mass
                md.db.dbr[i-1].num[j-1] = md.db.dbr[i-1].num[j-1]*(x0 - mrm)/x0;
                md.db.dbr[i-1].mass[j-1] = x0 - mrm;
                md.ml.mpp.enem = (md.ml.mpp.enem*md.ml.mass + erm*mrm)/(md.ml.mass + mrm);
                md.ml.mass = md.ml.mass + mrm;
                // ! mass of the debris stack node can become zero
                // !     => checked and collapsed after j loop is finished
                mrmtot = mrmtot + mrm;
            }

            // ! heat release => hex buffer, debris heat release total
            depomhex(md.hx, time, time, qrel); // ! deposit in heat exchange buffer
            md.db.qreltot = md.db.qreltot + qrel;

            // ! update debris bed stack node properties for the new temperature
            if(md.db.dbr[i-1].mass[j-1] > 0.0) {
                md.db.dbr[i-1].mpp[j-1].tm = tavn;
                md.db.dbr[i-1].qrel[j-1] = qrel/dt; // ! heat release stored as power (W)
                getmprop(md.mppb, md.db.dbr[i-1].mpp[j-1]);
                x0 = md.db.dbr[i-1].mass[j-1]/md.db.dbr[i-1].mpp[j-1].rhom; // ! debris volume in a node
                x1 = md.db.dbr[i-1].num[j-1]; // ! number of particles, preserved
                md.db.dbr[i-1].dia[j-1] = pow((6.0/pi*x0/x1),(1.0/3.0));
                // ! md.db.dbr[i-1].ars[j-1] = 6.0*x0/md.db.dbr[i-1].dia[j-1]
                md.db.dbr[i-1].ars[j-1] = pow((36.0*pi*pow(x0,2)*x1),(1.0/3.0));
                // ! update surface temperature (with old dia and props.)
                // ! for stabilization, qflx=h*(tsf-tsat), x2=0.1*h*dia/lam
                x2 = 0.1*qflx/max(tsf-tsat,1.0)*md.db.dbr[i-1].dia[j-1]/md.db.dbr[i-1].mpp[j-1].lamm;
                tsf = (tavn + x2*tsat)/(1.0 + x2);
                tsf = min(max(tsf,0.5*(tavn+tsat)),tavn);
                md.db.dbr[i-1].tsf[j-1] = tsf;
            }
        }

        // ! check the stack of ring i and erase mass=0 nodes
        // ! the number of nodes, md.db.dbr[i-1].ndz can change in the loop
        // j = ndz
        for (j=ndz;j==1;j--) // do
        {
            if(md.db.dbr[i-1].mass[j-1] < 1e-6) {dbstackcollapse(md.db.dbr[i-1],j);}
            // j = j - 1;
        }
        // if(j==0) exit 
        // end do
    }
    // << juwook 2022
    md.db.jj_g = md.db.jj_g/Adhf; 
    md.mp.qfluxmean = (md.db.qrel + md.ml.qrel)/Atot;
    // ! save the mass of debris remelted
    md.db.masstolmp = md.db.masstolmp + mrmtot;

    // ! adjust vertical boundaries of debris stack nodes
    fixdbednode(md.db);

    // ! check debris bed profile for avalanche
    if(md.db.irep > 0) {dbrepose(md.db,md.mppb);}

    // ! consider late in-side avalanche:
    // ! the avalanche into the in-side direction has a resistance by the
    // ! adjacent particles making a "ring" together. different from the case
    // ! of arriving particles, particles already in a debris bed may not 
    // ! move inside easily. thus, a separate switch and criterion is given
    // ! for the late in-side avalanche. 
    if(md.db.irepinlate > 0) {dbrepose_in(md.db,md.mppb);}

    // ! melt lump energy balance: phase 2
    // ! correct for heat release
    if(md.ml.mass > 0.0 && qbtmtot > 0.0){
        x0 = qbtmtot*dt; // ! total heat release (J)
        x1 = x0/md.ml.mass; // ! per mass
        md.ml.mpp.enem = md.ml.mpp.enem - x1;
        menetemp(md.mppb,md.ml.mpp.enem,md.ml.mpp.tm);
        // ! check undershoot (tm >= tl)
        if(md.ml.mpp.tm < md.clp.tl) {md.ml.mpp.tm = md.clp.tl;}
        // ! update melt phys. props. for new temperature
        getmprop(md.mppb,md.ml.mpp);

        // ! heat release => hex buffer, lump heat release total
        depomhex(md.hx, time, time, x0); // ! deposit in heat exchange buffer
        md.ml.qreltot = md.ml.qreltot + x0; // ! total heat release of lump (J)
        md.ml.qrel = qbtmtot; // ! heat release rate of lump (W)
    }

    // ! consider decay heat remaining in other places (e.g. the reactor vesse)
    // ! as a direct heat source in the containment => deposit to hex buffer 
    if( iresdh > 0 ){
        // !x0 = qdecaymass * (coremass - md.mj.massintot) *dt
        x0 = qdecaymass * (coremass/(1.0-fqdcex) - md.mj.massintot) *dt;
        depomhex(md.hx, time, time, x0); // ! deposit in heat exchange buffer
    }
}
/*###############################################################################*/ 




/*###############################################################################*/ 
void readmhex(melthex &mh, double time,double dt, double &qrate)
{
    int i0, i1, i;
    double q, qa; // ! buffer to count heat (J)
    
    q = 0.0; 
    qa = 0.0;

    // ! time is in thex(i-1)..thex(i)
    i0 = ceil(time/mh.dthex); // ! index for time
    i1 = ceil((time+dt)/mh.dthex); // ! for time+dt

    // ! do the 1st bin
    if(i0 > 0 && i0 <= nthexmax){
        q = mh.qhex[i0-1]; // ! rate (W) for 1st bin
        qa = mh.dqthestep*(time-mh.thex[i0-1]); // ! additional heat for the "current" bin (J)
                                            // ! =>see compensation... above
    }

    if(i1 == i0){
        q = q + qa / dt;
    }
    else {
        // ! need to cover multiple bins; do with integrated heat (W)
        q = q*(mh.thex[i0]-time);
        // ! do the last bin
        if(i1 > 0 && i1 <= nthexmax) {
            q = q + mh.qhex[i1-1]*(time+dt-mh.thex[i1-1]);
        }
        // ! do the intermediate bins
        for (i=i0+1; i<=i1-1; i++) // do i=i0+1, i1-1, 1
        {
            if(i > 0 && i <= nthexmax) {
                q = q + mh.qhex[i-1]*(mh.thex[i]-mh.thex[i-1]);
            }
        }
        q = (q+qa)/dt;
    }

    qrate = q;

    if(time > mh.timedq) {mh.dqthestep = 0.0;} // ! clear dqthestep when time passed
}
/*###############################################################################*/ 




/*###############################################################################*/ 
void getdbeddat (debrbed db, mpropbase mppb, double &mtot, double &vtot, double &etot, 
                 double &eav, double &tav, double &arstot, double &numtot, double &dsa, 
                 double &dmm, double &dmmlin, double &htop, double &hav, double &rad, 
                 double &tpk, double &rtpk, double &ztpk)
{
    double dmass[ndrndzmax], dmasso[ndrndzmax]; //, ddia[ndrndzmax]
    double hmass, m0, m1, d0, d1, dbvol;
    int i, j, nmax; //, ord[ndrndzmax];
    vector<double> ddia = vector<double>(ndrndzmax, 0.0);
    vector<int> ord = vector<int>(ndrndzmax, 0);

    // ! clear sums
    mtot=0.0; 
    vtot=0.0; 
    etot=0.0; 
    arstot=0.0; 
    numtot=0.0;

    // ! detect peak temperatur
    tpk = 0.0; 
    rtpk = 0.0; 
    ztpk = 0.0;

    // ! copy mass,dia data to work buffer, also calc. the sums
    nmax = 0;
    dbvol = 0.0;
    htop = 0.0;
    rad = 0.0;

    for (i=1; i<=db.ndr; i++) // do i = 1,db.ndr
    {
        if(db.dbr[i-1].mass[0] > 1e-3) {rad = db.rb[i];} // ! detect meaningful debris bed radius (>1mm thick)
        htop = max(htop, db.dbr[i-1].zb[db.dbr[i-1].ndz]);
        for (j=1; j<=db.dbr[i-1].ndz; j++) // do j = 1,db.dbr[i-1].ndz
        {
            if(db.dbr[i-1].mass[j-1] > 0.0){
                nmax = nmax + 1;
                dmass[nmax-1] = db.dbr[i-1].mass[j-1];
                ddia[nmax-1] = db.dbr[i-1].dia[j-1];
                ord[nmax-1] = nmax;
        
                mtot = mtot + db.dbr[i-1].mass[j-1];
                vtot = vtot + db.dbr[i-1].mass[j-1]/db.dbr[i-1].mpp[j-1].rhom;
                etot = etot + db.dbr[i-1].mass[j-1] * db.dbr[i-1].mpp[j-1].enem;
                arstot = arstot + db.dbr[i-1].ars[j-1];
                numtot = numtot + db.dbr[i-1].num[j-1];

                dbvol = dbvol + db.dar[i-1]*db.dbr[i-1].dz[j-1];

                // ! check peak temperature
                if(db.dbr[i-1].mpp[j-1].tm > tpk){
                    tpk = db.dbr[i-1].mpp[j-1].tm;
                    rtpk = db.rc[i-1];
                    ztpk = db.dbr[i-1].zc[j-1];
                }
            }
        }
    }

    if(rad > 0.0){
        hav = dbvol/(pi*pow(rad,2));
    }
    else {
        hav = 0.0;
    }

    if(mtot > 0.0){
        eav = etot/mtot;
    }
    else {
        eav = 0.0;
    }
    menetemp(mppb,eav,tav);
    etot = max(etot - mtot*enem0base, 0.0);
    if(arstot > 0.0){
        dsa = 6.0*vtot/arstot;
    }
    else {
        dsa = 0.0;
    }

    // ! calc. dmm
    // ! sort
    if(nmax > 10){
        qsort(ddia,1,1,nmax,ord);
    }
    else if(nmax > 1){
        // ssort(ddia,1,1,nmax,ord);
    }
    // ! re-order mass array as sorted, and find the half mass point
    hmass = mtot/2.0;
    m1 = 0.0;
    for (i=1; i<=nmax; i++) // do i = 1,nmax
    {
        dmasso[i-1] = dmass[ord[i-1]-1];
        m0 = m1;
        m1 = m1 + dmasso[i-1];
        if(hmass < m1) {break;}
    }
    if(i == 1){
        dmm = ddia[0]; 
        dmmlin = ddia[0];
    }
    else if(i > nmax){
        dmm = ddia[nmax-1]; 
        dmmlin = ddia[nmax-1];
    }
    else {
        d0 = ddia[i-2]; 
        d1 = ddia[i-1];
        dmm = d0 * pow((d1/d0),((hmass-m0)/(m1-m0)));
        dmmlin = d0 + (d1-d0)*(hmass-m0)/(m1-m0);
    }
}
/*###############################################################################*/ 













        



/*###############################################################################*/ 
// ! output dbrcool components data
// ! can be called only by the "mode" as i (initialize), w (write) or c (close)
void outdbrcool(string mode, string odir, meltdebr &md, double time=0,int step=0)//,ist)
{
    // character??

    double m,e,em,eml,emd,ehx,ta,ar,num,dsa,dmm,dmml,v,q,//qtothex,
        dbhtop,dbhav,dbrad,mlvol,tpk,rtpk,ztpk;
    int i, j, ichk, stp=1;//, ihexout, ioutlidx;

    static int ihexout;
    static double qtothex;
    static int ioutlidx;


    // = = = = = = = =
    // mode = "w";
    ofstream ojh, obh, olh, odh, ohh, odl, opr; // opt, 
    ojh.precision(4);
    obh.precision(4);
    olh.precision(4);
    odh.precision(4);
    ohh.precision(4);
    odl.precision(4);
    // opt.precision(4);
    opr.precision(4);
    // = = = = = = = =

    if (mode=="i"){ // 1646: mode = "i"
        
        printf("write - init. \n");
        system("mkdir mbrq");

        ojh.open("mbrq/ojh",ios::out);
            ojh << "# melt jet history; melt total internal energy basis=300K" << endl;
            ojh << "#1:time 2:d0 3:v0 4:temp0 5:di 6:vi 7:db 8:enem 9:hhf 10:hhp 11:zjl 12:llbr 13:llj 14:massin 15:enein" << endl;
            ojh << "\t" << fixed << time << "\t" << scientific << md.mj.d0 << "\t" << md.mj.v0<< "\t" 
                <<md.mj.temp0<< "\t" <<md.mj.di<< "\t" <<md.mj.vi<< "\t" 
                <<md.mj.db<< "\t" <<md.mj.mpp.enem<< "\t" <<md.plg.hhf<< "\t" 
                <<md.plg.hhp<< "\t" <<md.mj.zjl<< "\t" <<md.mj.llbr<< "\t" 
                <<md.mj.llj<< "\t" <<md.mj.massintot<< "\t" <<md.mj.eneintot<< endl;
        ojh.close();

        obh.open("mbrq/obh", ios::out);
            obh << "# mass/energy balance history; total mass (kg), cumulative energy(J); "
                   "melt total energy basis 300K; and coolant condition" << endl;
            obh << "#1:time 2:massmjetin 3:masspartot 4:massmlump 5:massdbed "
                   "6:enemjetin 7:enemlump 8:enedbed 9:hextot 10:qreltot_db 11:qreltot_lmp "
                   "12:qdctot_db 13:qdctot_lmp "
                   "14:wpoolmass 15:wpoolmassevap 16:wpoolenth 17:wpooleneevap 18:tl 19:tsat 20:pres "
                   "21:totqsrcex 22:massdbrtolmp 23:qrel_db 24:qrel_lmp 25:qcd_db 26:qdc_lmp" << endl;
            obh << "\t" << fixed << time << "\t" << scientific << md.mj.massintot<< "\t" << md.masspartot<< "\t" << md.ml.mass << "\t" 
                << "\t" << m << "\t" << md.mj.eneintot<< "\t" << eml<< "\t" << emd<< "\t" << ehx<< "\t" 
                << md.db.qreltot<< "\t" << md.ml.qreltot<< "\t" << md.db.qdctot<< "\t" << md.ml.qdctot<< "\t" 
                << md.wpoolmass<< "\t" <<md.wpoolmassevap<< "\t" <<md.wpoolenth<< "\t" <<md.wpooleneevap<< "\t" 
                << md.clp.tl<< "\t" <<md.clp.tspt<< "\t" <<md.clp.ptot<< "\t" <<md.totqsrcex<< "\t" <<md.db.masstolmp<< "\t" 
                << md.db.qrel<< "\t" <<md.ml.qrel<< "\t" <<md.db.qdc<< "\t" <<md.ml.qdc << endl;
            obh.close();
        obh.close();

        olh.open("mbrq/olh", ios::out);
            olh << "# melt lump history; melt total internal energy basis=300K; thickness is for heat transfer estimation (min.lim=1cm)\n" \
                   "#1:time 2:mass 3:ene 4:enetot 5:temp 6:volume 7:thickness 8:massfromdebris" << endl;
            olh << "\t" << fixed << time << "\t" << md.ml.mass << "\t" << e << "\t"
                << eml << "\t" << md.ml.mpp.tm << "\t" << mlvol << "\t" << md.ml.thk
                << "\t" << md.db.masstolmp << endl;
        olh.close();

        odh.open("mbrq/odh", ios::out);
            odh << "# debris bed history; melt total internal energy basis=300K\n" \
                   "#1:time 2:mass 3:ene 4:enetot 5:tempav 6:arstot 7:num \
                   8:diasau 9:diamm 10:diammlin 11:dmmbase 12:htop 13:have 14:radius \
                   15:tpeak 16:rtpeak 17:ztpeak" << endl;
            odh << "\t" << fixed << time << "\t" << scientific << m << "\t" << e << "\t" << emd << "\t" 
                << ta << "\t" << ar << "\t" << num << "\t" << dsa << "\t" << dmm << "\t" 
                << dmml << "\t" << md.mp.dmmb << "\t" << dbhtop << "\t" << dbhav << "\t" 
                << dbrad << "\t" << tpk << "\t" << rtpk << "\t" << ztpk << endl;
        odh.close();

        ohh.open("mbrq/ohh", ios::out);
            ohh << "# heat exchange history; qrate(W) q(J)(section) qtot(J)(cumulative)\n"\
                   "#1:id 2:tstart 3:tend 4:qrate 5:q 6:qtot" << endl;
        ohh.close();
        // initialization... whang
        ihexout = 1;
        qtothex=0.0;
        ioutlidx=0;


        opr.open("mbrq/opr", ios::out);
            opr << "# particle summary records"\
                   "#1:id 2:tgen 3:zgen 4:tset 5:rset 6:irset 7:mass 8:num \
                    9:diaset 10:eneset 11:tempset 12:tsfset" << endl;
        opr.close();



        odl.open("mbrq/odl", ios::out);
            odl << "" << endl;
        odl.close();

    }
    else if (mode=="w" || mode == "W"){ // 1684: mode = "w"
        // if not present?

        // ! write hex data when a section of time is passed
        if(time > md.hx.thex[ihexout]){
            do{
                q = md.hx.qhex[ihexout-1]*(md.hx.thex[ihexout]-md.hx.thex[ihexout-1]);
                qtothex = qtothex + q;
                ohh.open("mbrq/ohh", ios::app);
                ohh << "\t" << fixed << ihexout << "\t" << md.hx.thex[ihexout-1] << "\t" 
                    << md.hx.thex[ihexout] << "\t" <<  scientific << md.hx.qhex[ihexout-2] << "\t" 
                    << q << "\t" << qtothex << endl;
                ohh.close();
                ihexout = ihexout + 1;
            } while ((time > md.hx.thex[ihexout]));
        }

        // 1700 ! history data
        if ((time >=md.touth && md.dtouth > 0) || stp ==0 || mode == "W") {
            ojh.open("mbrq/ojh",ios::app);
            ojh << "\t" << fixed << time << "\t" << scientific << md.mj.d0 << "\t" << md.mj.v0<< "\t" 
                <<md.mj.temp0<< "\t" <<md.mj.di<< "\t" <<md.mj.vi<< "\t" 
                <<md.mj.db<< "\t" <<md.mj.mpp.enem<< "\t" <<md.plg.hhf<< "\t" 
                <<md.plg.hhp<< "\t" <<md.mj.zjl<< "\t" <<md.mj.llbr<< "\t" 
                <<md.mj.llj<< "\t" <<md.mj.massintot<< "\t" <<md.mj.eneintot<< endl;
            ojh.close();

            e = md.ml.mpp.enem;
            eml = md.ml.mass*max(e-enem0base, 0.0);
            if(md.ml.mass > 0.0){
                mlvol = md.ml.mass/md.ml.mpp.rhom;
                olh.open("mbrq/olh",ios::app); 
                olh << "\t" << fixed << time << "\t" << scientific << md.ml.mass << "\t" << e << "\t"
                    << eml << "\t" << md.ml.mpp.tm << "\t" << mlvol << "\t" << md.ml.thk
                    << "\t" << md.db.masstolmp << endl;
                olh.close();
            }

            getdbeddat(md.db,md.mppb, m,v,emd,e,ta,ar,num,dsa,dmm,dmml,
                       dbhtop,dbhav,dbrad,tpk,rtpk,ztpk);
            if(m>0.0){
                odh.open("mbrq/odh",ios::app);
                odh << "\t" << fixed << time << "\t" << m << "\t" << e << "\t" << emd << "\t" 
                    << ta << "\t" << ar << "\t" << num << "\t" << dsa << "\t" << dmm << "\t" 
                    << dmml << "\t" << md.mp.dmmb << "\t" << dbhtop << "\t" << dbhav << "\t" 
                    << dbrad << "\t" << tpk << "\t" << rtpk << "\t" << ztpk << endl;
                odh.close();
            }
            ehx = 0.0;
            for (i=1; i<=ihexout-1; i++){
                ehx = ehx + md.hx.qhex[i-1]*(md.hx.thex[i]-md.hx.thex[i-1]);
            }
            ehx = ehx + md.hx.qhex[ihexout-1]*(time-md.hx.thex[ihexout-1]); // ! hex till now

            obh.open("mbrq/obh",ios::app);
            obh << "\t" << fixed << time << "\t" << md.mj.massintot<< "\t" << md.masspartot<< "\t" << md.ml.mass << "\t" 
                << "\t" << m<< "\t" << md.mj.eneintot<< "\t" << eml<< "\t" << emd<< "\t" << ehx<< "\t" 
                << md.db.qreltot<< "\t" << md.ml.qreltot<< "\t" << md.db.qdctot<< "\t" << md.ml.qdctot<< "\t" 
                << md.wpoolmass<< "\t" <<md.wpoolmassevap<< "\t" <<md.wpoolenth<< "\t" <<md.wpooleneevap<< "\t" 
                << md.clp.tl<< "\t" <<md.clp.tspt<< "\t" <<md.clp.ptot<< "\t" <<md.totqsrcex<< "\t" <<md.db.masstolmp<< "\t" 
                << md.db.qrel<< "\t" <<md.ml.qrel<< "\t" <<md.db.qdc<< "\t" <<md.ml.qdc << endl;
            obh.close();

            md.touth = md.touth + md.dtouth;
        }

        // 1735 list data

        if(((time >= md.toutl && md.dtoutl > 0.0) || stp == 0) && m>0.0) {
            odl.open("mbrq/odl", ios::app);
            odl << "#index: " << ioutlidx << endl;
            odl << "#time= " << time << " step= " << stp << " (a blank-line delimits debris columns, two delimits records)" << endl;
            odl << "#1:i 2:iz 3:rc 4:zc 5:dr 6:dz 7:mass 8:ene 9:enetot 10:tempav 11:arstot 12:num 13:dia(sau) 14:tsf 15:qrel 16:ihtmod" << endl;
            for (i=1; i<=md.db.ndr; i++){
                ichk = 0;
                for (j=1; j<=md.db.dbr[i-1].ndz; j++){
                    e = md.db.dbr[i-1].mpp[j-1].enem;
                    em = md.db.dbr[i-1].mass[j-1] * max(e-enem0base, 0.0);
                    menetemp(md.mppb,e,ta);
                    odl << "\t" << fixed << i << "\t" << j << "\t" << scientific << md.db.rc[i-1] << "\t" << md.db.dbr[i-1].zc[j-1] << "\t" << md.db.rb[i]-md.db.rb[i-1] << "\t" 
                        << md.db.dbr[i-1].dz[j-1] << "\t" <<  md.db.dbr[i-1].mass[j-1] << "\t" << e << "\t" << em << "\t" << ta << "\t" 
                        << md.db.dbr[i-1].ars[j-1] << "\t" << md.db.dbr[i-1].num[j-1] << "\t" << md.db.dbr[i-1].dia[j-1] << "\t" << md.db.dbr[i-1].tsf[j-1] << "\t" << md.db.dbr[i-1].qrel[j-1] << "\t" 
                        << fixed << md.db.dbr[i-1].ihtmod[j-1] << endl;
                    ichk = 1;
                } // end do
                if(ichk == 1) {odl << "" << endl;}
            } // end do
            if(ichk == 1) {odl << "" << endl;}
            odl << "" << endl;

            md.toutl = md.toutl + md.dtoutl;
            ioutlidx = ioutlidx + 1;
        }
        // end if
    }
    else { // ! close mode "c"
        obh.close();
        odh.close();
        ohh.close();
        ojh.close();
        olh.close();

        odl.close();
    }
}
/*###############################################################################*/ 


/*###############################################################################*/ 
void outparsieve (string odir, meltdebr md)
{
    int u, i ;
    double cum=0.0, cumgen=0.0, mtotgen, mtot;

    ofstream ops;
    ops.precision(4);

    ops.open("mbrq/ops", ios::out);
    ops << "# particle sieve data (_gen:at generation, _solsed:for solid sediment)" << endl;
    ops << "#1:dia_high 2:mass_gen 3:massfrac_gen 4:cumfrac_gen 5:mass_solsed 6:massfrac_solsed 7:cumfrac_solsed" << endl;

    mtotgen = max(md.psievegen.masstot,1e-6);
    mtot = max(md.psieve.masstot,1e-6);

    for (i=1; i<=md.psieve.nn; i++){
        cumgen = cumgen + md.psievegen.mass[i-1]/mtotgen;
        cum = cum + md.psieve.mass[i-1]/mtot;
        ops << "\t" << scientific << md.psievegen.diab[i-1] << "\t" << md.psievegen.mass[i-1] << "\t" 
            <<  md.psievegen.mass[i-1]/mtotgen << "\t" <<  cumgen << "\t" << md.psieve.mass[i-1] << "\t" << 
            md.psieve.mass[i-1]/mtot << "\t" <<  cum << endl;
    } // end do

    ops << " " << endl;
    ops << "# total_mass_gen:    " << "\t" << md.psievegen.masstot << endl;
    ops << "# total_mass_solsed: " << "\t" << md.psieve.masstot << endl;

    ops.close();
}
/*###############################################################################*/ 



#endif