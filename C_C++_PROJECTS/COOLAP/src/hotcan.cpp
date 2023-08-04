#ifndef HOTCAN_H
#define HOTCAN_H

using namespace std;
#include <algorithm>
#include <vector>

#include "inttable_mod.cpp"
// #include "wrsttab.cpp"
#include "consts.cpp"
#include "idgas.cpp"
#include "math_mod.cpp"
/*###############################################################################*/



/*###############################################################################*/
const int nmaxqsrc=50, nmaxaelev=50, nnewloop=20; 
double pliml=614.0, plimh=1e8, tliml=274.0, tlimh=3000.0;

struct hhole
{
    double h=0.0, a=0.0, u=0.0, fa=1.0; 
    timetab tfa;
};

struct relvalve
{
    double a=0.0, po=0.0, pc=0.0, fa=0.0;
};

struct htcan
{
    int igiwrk=0 , // ! switch for gas inlet work evaluation
                // ! 0: with ptot; separate (injection work calculated with sum of (mass flow)/(density))
                // ! 1: with ptot; mixed (injection work calculated with (total mass flow)/(total density))
                // ! 2: with outside assumed pressures pi*; separate (rhoi*;ei* are evaluated at pi*;ti*)
        ifon=0,  // ! switch for consideration of non-uniformity of 
                // !  non-condensible gas "n" in the exit flow rate.
                // ! 0: not condider it
                // ! 1: consider it with time constant tau=fontau*(mdin/rhon)/volg and 
                // !    asymptotic ratio fonlim < 1
        nqsrc=0,  // ! number of the time dependent heat sources
        naelev=0, // ! number of the height dependent leak holes => aelev()
        ilossdp=0;  // ! switch for the method of loss consideration
                    // ! 0: give mass flow rate by mdloss* (or their time tables) (default)
                    // ! 1: give it by (ptot-penv)*mdloss* (pressure dependent flow out)
    double vol=0.0,   // ! total volume of the can (m3)
        volg=0.0, voll=0.0, vollc=0.0,   // ! volumes of gas, water, sub-cooled water (m3)
                // ! g: gas phase (vapor and non-condensibles) at Tsat(pv)
                // ! l: water at saturation temp.
                // ! lc: water staying subcooled at the bottom
        ts,tl,tlc,   // ! ts=tg,tn,tni=tsat(pv) (K), tl=tsat(pv'), pv'=pv-(fnr-1)*(pn+pni)
                // ! when fnr>1, water is in the equilibrium at pv'<pv (see fnr below)
        pv=0.0, pn=0.0, pni=0.0, ptot,    // ! pressures:vapor,non-condensibles, total (Pa)
                // ! v: vapor at equilibrium with water
                // ! n: non-condendible gas (can be non-uniformly mixed with vapor)
                // !        (n2 injected is assumed)
                // ! ni: non-condensible has (always uniformly mixed with vapor)
                // !        (internally generated gases, e.g. h2, are assumed)
        pout=1e5,    // ! pressure outside the holes (Pa) (default: atmospheric)
        poutlev=0.0, // ! level of the outlet connection IN THE NEXT CAN, for the calculation of
                        // ! static head for exhaust (should be initialized accordingly)
        penv=1e5,   // ! pressure of environment (outer world) (Pa) used to get the mass "loss"
                    // ! when ilossdp=1 is set
        fnr=1.0,    // ! non-condensible gas rich factor (pn'=fnr*pn pni'=fnr*pni 
                    // ! at water surface due to condensation; it shifts the equilibrium so that
                    // ! pv becomes higher, and models condensation retard by non-condensible gases)
        mv, mn, mni, ml, mlc,  // ! masses:vapor,non-condensibles,water,subcooled water(kg)
        mva, mla,  // ! intermediate values of the mass of vapor and water (after massevol)
        dmwev,   // ! evaporation(+)/condensation(-) mass change (kg) in a time step
        mmv=18e-3, mmn=28e-3, mmni=2e-3,  // ! molecular weight of the gases v, n, ni (kg/mol)
        ccv=4.45, ccn=3.5, ccni=3.5,  // ! cp*mm/rgasgen of the gases v, n, ni
                // ! cp*mm/rgasgen = 3.5 for two-atom molecules, about 4.45 for H2O at 1e5Pa, 380K
        rhov,rhon,rhoni,rhol,rholc, // ! densities:vapor,non-condensibles,water,subcooled water(kg/m3)
                // ! densities of gas components are at their partial pressures
        ev,en,eni,el,elc, // !internal energes:vapor,non-condensibles,water,subcooled water(J/kg)
        rhotot,  // ! total density of the gas (=rhov+rhon+rhoni)
        cvv,cvn,cvni,cpv,cpn,cpni, // ! specific heats of gases (J/kgK) needed for flow rate calc.
        cpl, betal,  // ! specific heat (K/kgK), volume expansion coeff.(1/K) of water 
                        // ! for heat transfer evaluation
        etot,  // ! total internal energy in the can (J)
        hlc=0.0,hl=0.0,htop=0.0, // ! heights (at the surface of subcooled water, saturated water; the top)
                        // ! (at present just evaluated if "tvollev" is given)  
        qsrc=0.0,  // ! heat source (in or out, in is positive) (W)
        qsrcex=0.0, // ! external heat source (W) as an interface with other modules
        tsrcex=0.0,  // ! temperature of the external heat source (K) (set the dry zone T if it exists) 
        wki=0.0, wko=0.0, // ! inlet/exit work (W) of the present step saved for 
                            // ! interfacing with other cans
        filc=0.0,  // ! fraction of injected water partitioned to subcooled pool
        fillost=0.0,  // ! fraction of injected water actually lost before injected
                    // ! mdilc=filc*(1-fillost)*mdil0, mdil=(1-filc)*(1-fillost)mdil0
                    // !  => give mdil0 and filc, fillost in the input 
        dtlc=0.0,  // ! dtlc/dt arbitrarily given (K/s)
        mdiv=0.0, mdin=0.0, mdini=0.0, mdil0=0.0, mdil, mdilc,  // ! mass in-flows (kg/s)
            // ! the in-flows are set 0 unless given by input (single values here or time tables below)
        mdlex=0.0,  // ! mass exchange from "lc" to "l" (negative is the opposit direction) (kg/s)
                        // ! (single value here or a table below)
        rhoiv=0.0,rhoin=0.0,rhoini=0.0,rhoil=0.0,  // ! inlet densities (kg/m3)
        eiv=0.0,ein=0.0,eini=0.0,eil=0.0,  // ! inlet internal energies (J/kg)
        tiv=0.0,tin=0.0,tini=0.0,til=0.0,  // ! inlet temperatures (K)
        piv=0.0,pin=0.0,pini=0.0, // ! inlet gas pressures (Pa) referred when igiwrk=2
            // ! INLET PROPERTIES: specify mdi*, ti*; ei* are evaluated by ptot and ti*
            // ! (in case water/vapor temp. is not consistent with EOS, fixed to near saturation state 
            // !  at the total pressure: a bit superheated vapor or a bit subcooled water)
        aeg=0.0,ael=0.0,aelc=0.0,  // ! areas of exit holes in gas, liquid and subcooled liquid (m2)
        mdov=0.0, mdon=0.0, mdoni=0.0, mdol=0.0, mdolc=0.0,  // ! mass out-flows (kg/s)
        mdorcic=0.0,  // ! out-flows by RCIC (kg/s)
        fon, fonlim= 1.0,  // ! factor on non-condensible gas "n" density at exit, its max. limit
        fov,  // ! factor on vapor and non-condensible gas mixing with it
        fontau = 1.0,  // ! factor on the time constant of fon application
        tifon=0.0,  // ! start time (s) of "n" injection (needed in fon evaluation)
        uog, uol, uolc, prc, // ! exit flow velocities ("uog" is not the exit velcotiy, it is
                        // ! (exit mass flow rate)/((hole area)*(upstream density))
                        // ! and ciritical flow index (prc=pr/pc,critical if<=1)
        fmdolim=0.5,  // ! factor of limit of flowout velocity considering gas pressure
                        // ! decrease and possible evaporation mass
        // ftqsrc[nmaxqsrc]={1.0},  // ! modification factor for tqsrc()
        mdlossg=0.0,mdlossv=0.0,mdlossn=0.0,mdlossni=0.0,mdlossl=0.0,mdlosslc=0.0, 
                        // ! arbitrary mass loss rate of fluids (leak to the environment) (kg/s)
                        // ! gas mass flow-out rate is distrbuted to components by densities.
        totmlossv=0.0,totmlossn=0.0,totmlossni=0.0,totmlossl=0.0,totmlosslc=0.0,toteloss=0.0, 
                        // ! cumulative loss of masses and energy by assumed leaks
        sbmg=0.0; // ! not zero: if a hole is submerged, the area is reduced to 1/10
        
        vector<double> ftqsrc = vector<double>(nmaxqsrc,1.0);

        timetab tqsrc[nmaxqsrc],  // ! heat sources => sum is given to qsrc, 
                                // ! actual number of the sources should be given as nqsrc
            tftqsrc[nmaxqsrc],// ! time dependent modification factor for tqsrc()=>set ftqsrc
            tmdiv, tmdin, tmdini, tmdil0,   // ! inlet mass flow rates: v,n,ni,l0=l+lc => mdi*
                // ! for l and lc, lc0 should be given and it is partitioned to l and lc by
                // !   mdil%x = (1-filc)*mdil0%x, mdilc%x = filc*mdil0%x
            tmdlex, tfilc, tfillost, tfnr, tdtlc,  // ! => mdlex, filc, fillost, fnr, dtlc
            tmdorcic, ttil, ttin,  // ! => mdorcic, til, tin
            taeg,tael,taelc,   // ! area of exit => ae*
            tmdlossg,tmdlossl,tmdlosslc; // ! mass loss rates => mdloss*
        inttab tvollev;  // ! table for the water volume => level relation
                        // ! u=volume (m3), x=level (m)
        hhole aelev[nmaxaelev];  // ! level dependent discharge holes
        relvalve relv; // ! a relief valve modelv
}; // end type htcan

// ! a variable set to carry in-out of mass(water,non-condensibles) and energy
struct meio
{
    double mwio,mnio,mniio,eio;
};
/*###############################################################################*/


// ! subroutines to handle the temporal evolution of a hotcan
// ! - initcan(hc, vol,vl,vlc,ts,tlc,pn,pni) : easily give initial equilibrium
// !     condition of the can with physical properties filled
// ! - stigas(p,t,mm,cc, rho,e,cv,cp,drdp,drdt) : calculate 
// !     non-condensible gas state vars (ideal gas eos)
// ! - massenevol : advance masses and total energy
// ! - getequil : get the thermal equilibrium at the new time step
// ! - preslevel : calculate static head at given height


/*###############################################################################*/
void initcan(struct htcan &hc)
{
    // ! * IF mm and cc of non-condensible gases are DIFFERENT FROM
    // !   THE DEFAULTS (mmn=28d-3,mmni=2d-3,ccn=3.5,ccni=4.45),
    // !   SPECIFY THEM BEFOR CALLING THIS.
    // ! * IF fnr>1 is wanted, SET IT IN ADVANCE.
    // ! * Either of ts or pv can be omitted, then calculate the saturation state.
    // ! elements read or calculated : ts,pv
    // ! elements read : vol, vl,vlc,tlc, pn, pni
    double ts,pv, vol, vl,vlc,tlc, pn, pni;
    int ir, i;
    double p,r,h,drp,drt,dhp,dht, tl;

    // ! vol and htop by table data
    if(hc.vol == 0.0 && hc.tvollev.n > 1) {hc.vol = hc.tvollev.u[hc.tvollev.n-1];}
    if(hc.htop == 0.0 && hc.tvollev.n > 1) {hc.htop = hc.tvollev.x[hc.tvollev.n-1];}
    // ! non-zero water levels override vol* settings  // !km170531
    if( hc.hl+hc.hlc != 0.0 && hc.tvollev.n > 1)
    {
        itabevalr(hc.tvollev, hc.hlc, hc.vollc);
        itabevalr(hc.tvollev, hc.hl, hc.voll);
        hc.voll = hc.voll - hc.vollc;
    } 

    // ! evaluate and check vol*
    if(hc.vol == 0.0) {hc.vol = hc.volg+hc.voll+hc.vollc;}
    if(hc.volg == 0.0) {hc.volg = hc.vol-hc.voll-hc.vollc;}
    if(hc.voll == 0.0) {hc.voll = hc.vol-hc.volg-hc.vollc;}
    if(hc.vollc == 0.0) {hc.vollc = max(hc.vol-hc.volg-hc.voll,0.0);}
    if(hc.vol<0.0 || hc.volg<0.0 || hc.voll<0.0 || hc.vollc<0.0) 
    {
        printf("TMP-initcan\n");
        // write(*,'(a)') "initcan: invalid volume spec."
        // write(*,'(a,1p,4(1x,e12.4))') "  vol,volg,voll,vollc=",hc.vol,hc.volg,hc.voll,hc.vollc
    }

    vol = hc.vol; 
    vl = hc.voll; 
    vlc = hc.vollc;
    hc.volg = vol - (vl + vlc);
    if(hc.volg < 0.0) {
        printf("TMP-initcan\n"); // call errorhdr("initcan: Volume in-consistency Vg<0","cont")
    }

    // ! some valid temp should be given to "lc"
    tlc = hc.tlc;
    if(hc.tlc < tliml) {hc.tlc = tliml;}

    // ! water state vars
    // ! saturation
    // !  call steamv(ts,p,v,h,s,cp,cv,hrho,cc,beta,akap,gam,1,i,1) // ! vapor, sat@t
    // !  hc.rhov = 1.0/v; hc.ev = h-p*v
    // !  hc.cvv = cv; hc.cpv = cp; hc.ccv = cp*hc.mmv/rgasgen
    if(hc.pv <= 0.0) {hc.pv = psatwaterf(hc.ts);}
    if(hc.ts <= 0.0) {hc.ts = tsatwaterf(hc.pv);}
    pv = hc.pv; 
    ts = hc.ts;
    /* ir = */ wrsteamtab(pv,ts,1,r,drp,drt,h,dhp,dht,ir); // ! sat. vapor at ps,ts
    hc.rhov = r; 
    hc.ev = h - pv/r;
    hc.cpv = dht; 
    hc.cvv = dht-ts/pow(r,2)*pow(drt,2)/drp; 
    hc.ccv = dht*hc.mmv/rgasgen;

    // ! consideration of condensation retard by non-condensible gases
    // ! if fnr=1 no effect
    hc.fnr = max(hc.fnr, 1.0);
    if((hc.fnr-1.0)*(pn+pni) > pv-pliml) {hc.fnr = 1.0+(pv-pliml)/(pn+pni);}
    p = pv - (hc.fnr-1.0)*(pn + pni);
    // !  call steamv(tl,p,v,h,s,cp,cv,hrho,cc,beta,akap,gam,4,i,1) // ! water, sat@p
    // !  hc.rhol = 1.0/v; hc.el = h-p*v; hc.cpl = cp; hc.betal = beta
    tl = tsatwaterf(p);
    /* ir = */ wrsteamtab(p,tl,0, r,drp,drt,h,dhp,dht,ir); // ! sat. water at p
    hc.tl = tl;
    hc.rhol = r; 
    hc.el = h - p/r; 
    hc.cpl = dht; 
    hc.betal = -drt/r;
    // ! subcooled water
    // !  p = psat(tlc)
    // !  call steamv(tlc,p,v,h,s,cp,cv,hrho,cc,beta,akap,gam,0,i,0)
    // !  hc.rholc = 1.0/v; hc.elc = h-p*v
    p = psatwaterf(tlc);
    /* ir = */ wrsteamtab(p,tlc,0, r,drp,drt,h,dhp,dht,ir); // ! subc. water
    hc.rholc = r; 
    hc.elc = h - p/r;

    // ! noncondensible gas state vars
    pn = hc.pn; 
    pni = hc.pni;
    stigas1(pn,ts,hc.mmn,hc.ccn,hc.rhon,hc.en,hc.cvn,hc.cpn);
    stigas1(pni,ts,hc.mmni,hc.ccni,hc.rhoni,hc.eni,hc.cvni,hc.cpni);

    // ! totals
    hc.mv = hc.volg*hc.rhov;
    hc.ml = hc.voll*hc.rhol;  
    hc.mlc = hc.vollc*hc.rholc;
    hc.mn = hc.volg*hc.rhon;  
    hc.mni = hc.volg*hc.rhoni;
    hc.ptot = hc.pv + hc.pn + hc.pni;
    hc.rhotot = hc.rhov + hc.rhon + hc.rhoni;
    hc.etot = hc.volg*(hc.rhov*hc.ev + hc.rhon*hc.en + hc.rhoni*hc.eni) 
        + hc.voll*hc.rhol*hc.el + hc.vollc*hc.rholc*hc.elc;

    // ! check limits of numbers
    if(hc.nqsrc > nmaxqsrc) { // call errorhdr("initcan: nqsrc > nmaxqsrc; abort.","stop")
        printf("TMP-hotcan-initcan\n");
        exit(EXIT_FAILURE);
    }
    if(hc.naelev > nmaxaelev) { // call errorhdr("initcan: naelev > nmaxaelev; abort.","stop")
        printf("TMP-hotcan-initcan\n");
        exit(EXIT_FAILURE);
    }

    // ! check of table inputs
    if(hc.tmdiv.n > 1)  {ttabcheck(hc.tmdiv);}
    if(hc.tmdin.n > 1)  {ttabcheck(hc.tmdin);}
    if(hc.tmdini.n > 1)  {ttabcheck(hc.tmdini);}
    if(hc.tmdil0.n > 1)  {ttabcheck(hc.tmdil0);}
    if(hc.tmdlex.n > 1)  {ttabcheck(hc.tmdlex);}
    if(hc.tfilc.n > 1)  {ttabcheck(hc.tfilc);}
    if(hc.tfillost.n > 1)  {ttabcheck(hc.tfillost);}
    if(hc.tfnr.n > 1)  {ttabcheck(hc.tfnr);}
    if(hc.taeg.n > 1)  {ttabcheck(hc.taeg);}
    if(hc.tael.n > 1)  {ttabcheck(hc.tael);}
    if(hc.taelc.n > 1)  {ttabcheck(hc.taelc);}
    if(hc.tdtlc.n > 1)  {ttabcheck(hc.tdtlc);}
    if(hc.tmdlossg.n > 1)  {ttabcheck(hc.tmdlossg);}
    if(hc.tmdlossl.n > 1)  {ttabcheck(hc.tmdlossl);}
    if(hc.tmdlosslc.n > 1)  {ttabcheck(hc.tmdlosslc);}
    if(hc.tmdorcic.n > 1)  {ttabcheck(hc.tmdorcic);}
    if(hc.ttil.n > 1)  {ttabcheck(hc.ttil);}
    if(hc.ttin.n > 1)  {ttabcheck(hc.ttin);}
    for (i=1; i<=hc.nqsrc;i++) // do i = 1, hc.nqsrc
    {
        if(hc.tqsrc[i-1].n > 0)  {ttabcheck(hc.tqsrc[i-1]);}
        if(hc.tftqsrc[i-1].n > 0)  {ttabcheck(hc.tftqsrc[i-1]);}
    }
    for (i=1; i<=hc.naelev; i++) // do i = 1, hc.naelev
    {
        if(hc.aelev[i-1].tfa.n > 0)  {ttabcheck(hc.aelev[i-1].tfa);}
    }
    if(hc.tvollev.n > 1)
    {
        itabcheck(hc.tvollev);
        if(hc.tvollev.typ != 0)
        {
            printf("'# initcan: The volume-level table has a wrong type.'\n");
            printf("'#   type 0 (linear interpolation) is forced.'\n");
            hc.tvollev.typ = 0;
        }
    }
    if(hc.naelev > 0 && hc.tvollev.n == 0)
    {
        printf ("'initcan: The height dependent leak holes are difined without'\n");
        printf ("'    difining volume-level table.'\n");
    }

    // ! do the water level
    if(hc.tvollev.n>1) 
    {
        itabeval(hc.tvollev, hc.vollc, hc.hlc);
        itabeval(hc.tvollev, hc.voll+hc.vollc, hc.hl);
        itabeval(hc.tvollev, hc.vol, hc.htop);
    }

}
/*###############################################################################*/

/*###############################################################################*/
void massenevol(htcan &hc, double time, double dt, meio &xio) // int &ist,
{
    int i, ir;
    double q0, q1, pr, pc, gam, vv,rv,rn,rni, rt, 
           mi = 1e-12, tau, wrki, wrko, t,r,h,drp,drt,dhp,dht,cv,cp, 
           mwio,mnio,mniio,eio, aleg, dtsp,dhfg,uogmax, 
           vdlossg,mlossv,mlossn,mlossni,mlossl,mlosslc,eloss,dpenv;

    // if(present(ist)) ist = 0

    // ! mass and energy change other than calculated here is possible
    // ! (e.g. some external route) 
    // ! they are GIVEN BEFORE CALLING this.

    // ! do the time tabled conditions: heat soruces, inlet flows, exit areas
    // ! and external heat input
    // ! a single constant value heat source "qsrc" is overridden by time table
    // ! heat source(s) and external heat source 
    if(hc.nqsrc > 0){
        q0 = 0.0;
        for (i=1; i<=hc.nqsrc; i++) // do i = 1, hc.nqsrc
        {
            q1 = 0.0;
            if(hc.tqsrc[i-1].n > 0) {ttabeval(hc.tqsrc[i-1], time, q1);}
            if(hc.tftqsrc[i-1].n > 0) {ttabeval(hc.tftqsrc[i-1], time, hc.ftqsrc[i-1]);}
            q0 = q0 + q1*hc.ftqsrc[i-1];  // ! tqsrc can be modified by given factors
        }
        hc.qsrc = q0;
        hc.qsrc = hc.qsrc + hc.qsrcex; // ! add external heat source if given
    }
    else if (hc.qsrcex != 0.0){
        // ! external heat source overrides the single qsrc value
        // ! qsrc is not updated if given directly (not by time table),{ adding 
        // ! qsrcex would make over counting at every time step.
        hc.qsrc = hc.qsrcex;
    }
    if(hc.tmdiv.n > 0) {ttabeval(hc.tmdiv, time, hc.mdiv);}
    if(hc.tmdin.n > 0) {ttabeval(hc.tmdin, time, hc.mdin);}
    if(hc.tmdini.n > 0) {ttabeval(hc.tmdini, time, hc.mdini);}
    if(hc.tmdil0.n > 0) {ttabeval(hc.tmdil0, time, hc.mdil0);}
    if(hc.tmdlex.n > 0) {ttabeval(hc.tmdlex, time, hc.mdlex);}
    if(hc.ml <= 0.0) {hc.mdlex = max(hc.mdlex, 0.0);}
    if(hc.mlc <= 0.0) {hc.mdlex = min(hc.mdlex, 0.0);}
    if(hc.tfilc.n > 0) {ttabeval(hc.tfilc, time, hc.filc);}
    if(hc.tfillost.n > 0) {ttabeval(hc.tfillost, time, hc.fillost);}
    if(hc.tfnr.n > 0) {ttabeval(hc.tfnr, time, hc.fnr);}
    hc.mdil = (1.0 - hc.filc)*(1.0 - hc.fillost)*hc.mdil0;
    hc.mdilc = hc.filc*(1.0 - hc.fillost)*hc.mdil0;
    if(hc.taeg.n > 0) {ttabeval(hc.taeg, time, hc.aeg);}
    if(hc.tael.n > 0) {ttabeval(hc.tael, time, hc.ael);}
    if(hc.taelc.n > 0) {ttabeval(hc.taelc, time, hc.aelc);}
    if(hc.tdtlc.n > 0) {ttabeval(hc.tdtlc, time, hc.dtlc);}
    if(hc.tmdlossg.n > 0) {ttabeval(hc.tmdlossg, time, hc.mdlossg);}
    if(hc.tmdlossl.n > 0) {ttabeval(hc.tmdlossl, time, hc.mdlossl);}
    if(hc.tmdlosslc.n > 0) {ttabeval(hc.tmdlosslc, time, hc.mdlosslc);}
    if(hc.tmdorcic.n > 0) {ttabeval(hc.tmdorcic, time, hc.mdorcic);} // ! RCIC out-flow
    if(hc.ttil.n > 0) {ttabeval(hc.ttil, time, hc.til);} // ! injection water temp.
    if(hc.ttin.n > 0) {ttabeval(hc.ttin, time, hc.tin);} // ! injection nitrogen temp.
    for(i=1;i<=hc.naelev;i++) // do i = 1, hc.naelev
    {
        if(hc.aelev[i-1].tfa.n > 0) {ttabeval(hc.aelev[i-1].tfa, time, hc.aelev[i-1].fa);}
    }
    // ! relief valve model if exists
    if(hc.relv.a > 0.0){
        if(hc.ptot > hc.relv.po+hc.pout) {hc.relv.fa = 1.0;}
        if(hc.ptot < hc.relv.pc+hc.pout) {hc.relv.fa = 0.0;}
    }

    // ! calc. inlet material properties
    if(hc.igiwrk != 2){ // ! if not giving the outside pressure for inlet work eval.
        hc.piv = hc.ptot; 
        hc.pin = hc.ptot; 
        hc.pini = hc.ptot;
    }
    if(hc.mdiv > 0.0 && hc.piv > 0.0){
        // !     steamz(t,hc.piv,v,h,s, 2) // ! sat.steam at ptot
        t = tsatwaterf(hc.piv);
        if(hc.tiv < t) {hc.tiv = t+1.0;}  // ! a bit superheated condition if temp. is wrong
        // !     steamz(hc.tiv,hc.piv,v,h,s, 0) // ! superheated vapor
        // !     hc.rhoiv = 1.0/v; hc.eiv = h - hc.ptot*v
        wrsteamtab(hc.piv,hc.tiv,1,r,drp,drt,h,dhp,dht,ir); // ! super heated vapor
        hc.rhoiv = r; 
        hc.eiv = h - hc.piv/r;
    }
    if(hc.mdin > 0.0 && hc.pin > 0.0) {
        stigas1(hc.pin,hc.tin,hc.mmn,hc.ccn,hc.rhoin,hc.ein,cv,cp);
    }
    if(hc.mdini > 0.0 && hc.pini > 0.0) {
        stigas1(hc.pini,hc.tini,hc.mmni,hc.ccni,hc.rhoini,hc.eini,cv,cp);
    }

    if(hc.mdil > 0.0 || hc.mdilc > 0.0){
        // !     steamz(t,hc.ptot,v,h,s, 4) // ! sat.water at ptot
        t = tsatwaterf(hc.ptot);
        if(hc.til > t) {hc.til = t-1.0;}  // ! a bit subcooled condition
        // !     steamz(hc.til,hc.ptot,v,h,s, 0) // ! subcooled water
        // !     hc.rhoil = 1.0/v; hc.eil = h - hc.ptot*v
        wrsteamtab(hc.ptot,hc.til,0,r,drp,drt,h,dhp,dht,ir); // ! subcooled water
        hc.rhoil = r; 
        hc.eil = h - hc.ptot/r;
    }

    // ! calc. exit flow rates
    // ! gas (non-uniformity of "n" is considered)
    aleg = 0.0; // ! area of gas phase leak hole specified by the level
    // ! relief valve model
    if(hc.relv.a > 0.0){
        aleg = aleg + hc.relv.a*hc.relv.fa;
    }
    // ! height dependent holes
    for(i=1; i<=hc.naelev; i++) // do i=1, hc.naelev
    {
        if(hc.aelev[i-1].h > hc.hl) aleg = aleg + hc.aelev[i-1].a*hc.aelev[i-1].fa;
    }
    // ! evaluate gas flow out when hole exists and ptot>pout
    if(hc.aeg+aleg > 0.0 && hc.ptot > hc.pout){
        if(hc.ifon == 1){ // ! do the modification of gas densities at the exit
            tau = hc.volg / max(hc.mdin/max(hc.rhoin,mi), mi) * hc.fontau;
            tau = max(tau, 1.0); // ! arbitrary limitation, though..
            hc.fon = hc.fonlim*(1.0 - exp(-max(time-hc.tifon, 0.0)/tau));
            hc.fov = (hc.rhov/hc.mmv + hc.rhoni/hc.mmni + (1.0-hc.fon)*hc.rhon/hc.mmn) 
                /max(hc.rhov/hc.mmv + hc.rhoni/hc.mmni, mi);
            rv = hc.fov*hc.rhov; rn = hc.fon*hc.rhon; rni = hc.fov*hc.rhoni;
            // ! assume the gas "n" is thinner at the exit than the average density in the can
            // ! due to the diffusion and convection mixing time constant "tau"
        }
        else { // ! use the balk densities as they are
            rv = hc.rhov; 
            rn = hc.rhon; 
            rni = hc.rhoni;
        }
        rt = max(rv + rn + rni, mi);
        // !// ! the modified gas densities are used in the exit work below, too.
        pr = hc.pout/hc.ptot;
        gam = (hc.cpv*rv+hc.cpn*rn+hc.cpni*rni)/max(hc.cvv*rv+hc.cvn*rn+hc.cvni*rni, mi);
        pc = pow((2.0/(gam+1.0)),(gam/(gam-1.0)));
        hc.prc = pr/pc;
        if(pr > pc){ // ! sub-critical flow
            hc.uog = sqrt( 2.0*gam/(gam-1.0)*hc.ptot/rt*( pow(pr,(2.0/gam))-pow(pr,((gam+1.0)/gam))));
        }
        else { // ! critical flow (not dependent on pr)
            hc.uog = sqrt( hc.ptot/rt*gam*pow((2.0/(gam+1.0)),((gam+1.0)/(gam-1.0))));
        }
        // ! check excess pressure over the back pressure and limit the velocity
        // ! to prevent overshooting
        // ! ug*a*dt/volg*ptot = dp <= ptot-pout ... flow out of gas
        // ! -> ug < max(0,ptot-pout)*volg/(ptot*a*dt)
        // ! additional flow out by evaporation
        // !  (dtsat/dp)*dpv*cpl*voll*rhol/dhfg/rhov (m3)
        uogmax = (hc.ptot - hc.pout)/hc.ptot*hc.volg/(hc.aeg+aleg)/dt; // ! gas phase flow out
        // ! additional flow out due to evaporation 
        if(hc.voll > 0.0){
            dtsp = dtsdpwaterf(hc.pv);
            dhfg = hfgwaterf(hc.pv);
            uogmax = uogmax + 
                dtsp*(hc.ptot - hc.pout)*hc.pv/hc.ptot*hc.cpl*hc.voll*hc.rhol 
                /dhfg/hc.rhov/(hc.aeg+aleg)/dt;
        }
        hc.uog = min(hc.uog, hc.fmdolim*uogmax);
        // !{ do it 
        hc.mdov  = (hc.aeg+aleg)*rv* hc.uog;
        hc.mdon  = (hc.aeg+aleg)*rn* hc.uog;
        hc.mdoni = (hc.aeg+aleg)*rni*hc.uog;
        // ! uog is not actually the exit velocity. 
        // ! uog=(mass flow rate)/((exit area)*(upstream density))
        // ! the exit velocity is actually
        // ! uex = F/(A*rhotot')=mdov/(A*rhov')=mdon/(A*rhon')=mdoni/(A*rhoni')
        // ! where rho*' is the densities at the exit. because the composition does not
        // ! change between the exit and the upstream, the density ratios are the same,
        // ! rhov':rhon':rhoni'=rhov:rhon:rhoni. thus, the following also holds,
        // ! uog = F/(A*rhotot)=mdov/(A*rhov)=mdon/(A*rhon)=mdoni/(A*rhoni)
        // ! where uog is a "converted exit velocity" based on the upstream densities.
        for(i=1;i<=hc.naelev;i++) // do i=1, hc.naelev
        {
            if(hc.aelev[i-1].h > hc.hl) hc.aelev[i-1].u = hc.uog;
        }
    }
    else {
        hc.mdov = 0.0;
        hc.mdon = 0.0;
        hc.mdoni = 0.0;
        hc.uog = 0.0;
        hc.prc = 0.0;
    }
    // ! check the lower limit steam pressure
    // ! solution fails when pv<pliml; non-condensible rich gas flow out may cause this.
    // ! in that case, steam flow is blocked as a remedy; deformation of the physics 
    // ! by this treatment is not significant.
    if(hc.pv <= pliml*2.0) {hc.mdov = 0.0;}

    // ! liquids
    if(hc.ael > 0.0){
        hc.uol = sqrt(2.0*max(hc.ptot - hc.pout,0.0)/hc.rhol);
        hc.mdol = hc.ael*hc.rhol*hc.uol;
    }
    else {
        hc.uol = 0.0;
        hc.mdol = 0.0;
    }
    for (i=1; i<=hc.naelev; i++) // do i=1, hc.naelev
    {
        if(hc.aelev[i-1].h > hc.hlc && hc.aelev[i-1].h <= hc.hl && hc.ml > 0.0){
            hc.aelev[i-1].u = 
                sqrt(2.0*max(hc.rhol*grav*(hc.hl-hc.aelev[i-1].h)+hc.ptot-hc.pout,0.0)/hc.rhol);
            if(hc.sbmg != 0.0){
                hc.mdol = hc.mdol + hc.sbmg*hc.aelev[i-1].a*hc.aelev[i-1].fa*hc.rhol*hc.aelev[i-1].u;
            }
            else {
                hc.mdol = hc.mdol + hc.aelev[i-1].a*hc.aelev[i-1].fa*hc.rhol*hc.aelev[i-1].u;
            }
        }
    }

    if(hc.aelc > 0.0){
        hc.uolc = sqrt(2.0*max(hc.ptot - hc.pout,0.0)/hc.rholc);
        hc.mdolc = hc.aelc*hc.rholc*hc.uolc;
    }
    else {
        hc.uolc = 0.0;
        hc.mdolc = 0.0;
    }
    for (i=1; i<=hc.naelev; i++) // do i=1, hc.naelev
    {
        if(hc.aelev[i-1].h <= hc.hlc && hc.mlc > 0.0){
            hc.aelev[i-1].u = sqrt(2.0*max(hc.rhol*grav*(hc.hl-hc.hlc)+ 
                hc.rholc*grav*(hc.hlc-hc.aelev[i-1].h)+hc.ptot-hc.pout,0.0)/hc.rholc);
            hc.mdolc = hc.mdolc + hc.aelev[i-1].a*hc.aelev[i-1].fa*hc.rholc*hc.aelev[i-1].u;
        }
    }

    // ! mass evolution
    mwio = dt*(hc.mdiv - hc.mdov)+dt*(hc.mdil - hc.mdol - hc.mdorcic)+dt*(hc.mdilc - hc.mdolc);
    mnio = dt*(hc.mdin - hc.mdon);
    mniio = dt*(hc.mdini - hc.mdoni);
    hc.mva = hc.mv + dt*(hc.mdiv - hc.mdov); // ! intermediate values, dmev is reserved
    hc.mla = hc.ml + dt*(hc.mdil - hc.mdol - hc.mdorcic + hc.mdlex); // !  same as above
    hc.mlc = hc.mlc + dt*(hc.mdilc - hc.mdolc - hc.mdlex);
    hc.mn = hc.mn + mnio;
    hc.mni = hc.mni + mniio;

    // ! energy evolution
    // ! energy flow-in/out and work from/to outside other than formulated here 
    // ! can be given by added to the internal energies in advance.
    // ! to pass the work ocurr in the present step, wrki/wrko are SAVED AT RETURN.

    // ! arbitrarily given change of subcooled water temperature, limited at tl
    // ! needed here to refer the time related variables, tlc is just passed to getequil 
    hc.tlc = min(hc.tlc + hc.dtlc*dt, hc.tl);

    // ! outlet work (gas exit flow rates and densities above are referred) (W)
    wrko =  hc.ptot*( (hc.mdov+hc.mdon+hc.mdoni)/max(rt,mi) 
            + hc.mdol/hc.rhol + hc.mdorcic/hc.rhol + hc.mdolc/hc.rholc );

    // ! inlet work
    // ! densities at inlet given as input is at total pressure and inlet temperatures.
    if(hc.igiwrk == 1){ // ! mixture injection
        vv = hc.mdiv/max(hc.rhoiv,mi) + hc.mdin/max(hc.rhoin,mi) + hc.mdini/max(hc.rhoini,mi);
        vv = max(vv,mi);
        rv = hc.mdiv/vv; rn = hc.mdin/vv; rni = hc.mdini/vv;
        rt = max(rv + rn + rni, mi);
        wrki = hc.ptot*(hc.mdiv+hc.mdin+hc.mdini)/rt;
    }
    else { // ! separate injection; 
        // ! if igiwrk==2 piv,pin,pini are given at outside and corresponding rho,e are used
        wrki = hc.piv*hc.mdiv/max(hc.rhoiv,mi) + hc.pin*hc.mdin/max(hc.rhoin,mi) 
            + hc.pini*hc.mdini/max(hc.rhoini,mi);
    }
    if( hc.mdil > 0.0 || hc.mdilc > 0.0) {
        wrki = wrki + hc.ptot*( hc.mdil/max(hc.rhoil,mi) + hc.mdilc/max(hc.rhoil,mi) );
    }
    // ! energy evolution
    eio = dt*( wrki - wrko 
        + hc.mdiv*hc.eiv + hc.mdin*hc.ein + hc.mdini*hc.eini + (hc.mdil + hc.mdilc)*hc.eil 
        -(hc.mdov*hc.ev + hc.mdon*hc.en + hc.mdoni*hc.eni + (hc.mdol + hc.mdorcic)*hc.el + hc.mdolc*hc.elc) 
        + hc.qsrc );
    hc.etot = hc.etot + eio;

    // if(present(xio)){
    xio.mwio = mwio; 
    xio.mnio = mnio; 
    xio.mniio = mniio; 
    xio.eio = eio;
    // end if

    // ! save works
    hc.wki = wrki; 
    hc.wko = wrko; // ! inlet and exit works (W)

    // ! arbitrary mass and energy (including work) loss
    // ! every gas component 
    rt = hc.rhov + hc.rhon * hc.rhoni;
    vdlossg = hc.mdlossg / max(rt,mi);
    hc.mdlossv = vdlossg*hc.rhov;
    hc.mdlossn = vdlossg*hc.rhon;
    hc.mdlossni = vdlossg*hc.rhoni;
    mlossv = dt*hc.mdlossv;
    mlossn = dt*hc.mdlossn;
    mlossni = dt*hc.mdlossni;
    mlossl = dt*hc.mdlossl;
    mlosslc = dt*hc.mdlosslc;
    // ! take mloss* as proportionality factor of flow rate and dp
    if(hc.ilossdp == 1){
        dpenv = max(hc.ptot - hc.penv, 0.0);
        mlossv = mlossv * dpenv;
        mlossn = mlossn * dpenv;
        mlossni = mlossni * dpenv;
        mlossl = mlossl * dpenv;
        mlosslc = mlosslc * dpenv;
    }
    if(hc.pv < pliml*2.0 || hc.mva < mlossv) {mlossv = 0.0;}
    if(hc.mn < mlossn) {mlossn = 0.0;}
    if(hc.mni < mlossni) {mlossni = 0.0;}
    if(hc.mla < mlossl) {mlossl = 0.0;}
    if(hc.mlc < mlosslc) {mlosslc = 0.0;}
    eloss = mlossv*hc.ev + mlossn*hc.en + mlossni*hc.eni + mlossl*hc.el + mlosslc*hc.elc 
            + hc.ptot*(vdlossg + mlossl/hc.rhol + mlosslc/hc.rholc);

    // ! do the changes
    hc.mva = hc.mva - mlossv;
    hc.mla = hc.mla - mlossl;
    hc.mn = hc.mn - mlossn;
    hc.mni = hc.mni - mlossni;
    hc.mlc = hc.mlc - mlosslc;
    hc.etot = hc.etot - eloss;

    // ! save losses
    hc.totmlossv = hc.totmlossv + mlossv;
    hc.totmlossn = hc.totmlossn + mlossn;
    hc.totmlossni = hc.totmlossni + mlossni;
    hc.totmlossl = hc.totmlossl + mlossl;
    hc.totmlosslc = hc.totmlosslc + mlosslc;
    hc.toteloss = hc.toteloss + eloss;
}
/*###############################################################################*/


/*###############################################################################*/
void getequil(htcan &hc, int itlcst, int &ist)
{
    int itl, is, neps, inew, ir;
    double mi = 1e-12, epsnew=1e-6, res, 
    f1, f2, f3, f4, f5, f6, f10, f20, f30, f40, f50, f60, 
    rv0, ev0, rl0, el0, rg0, v0, t0=273.15, 
    dmev, ddmev, dpv, dpn, dpni, dvg, dvl, 
    aevpv, aelpv, aenpv, aenipv, arvpv, arlpv, 
    arnpv, arnpn, arnipv, arnipni, aelpn, arlpn, 
    ts, dtsp, tl, p,r,drp,drt,h,dhp,dht,  
    ts1, dts, aevt, aent, aenit, arvt, arnt, arnit, 
    pvorg,pnorg,pniorg,tsorg, ptotorg;

    // a(6,7), a1(4,5), 
    vector<vector<long double>> a = vector<vector<long double>>(6, vector<long double>(7,0.0));
    vector<vector<long double>> a1 = vector<vector<long double>>(4, vector<long double>(5,0.0));
    
    itl = 0;
    is = 0;

    itl = itlcst;  // ! if itlcst=1, tl must be set in advance.

    // ! backup the states for the case of running out of water
    pvorg = hc.pv; 
    pnorg = hc.pn; 
    pniorg = hc.pni; 
    tsorg = hc.ts; 
    ptotorg = hc.ptot;

    // ! default scale of the newton eqs
    // !  steamv(400.0,1e5,v,h,s,cp,cv,hrho,cc,beta,akap,gam,0,i,0) // ! standard vapor
    // !  rv0 = 1.0/v; ev0 = h
    // !  steamv(300.0,1e5,v,h,s,cp,cv,hrho,cc,beta,akap,gam,0,i,0) // ! standard water
    // !  rl0 = 1.0/v; el0 = h-v*1e5
    wrsteamtab(1e5,400.0,1,r,drp,drt,h,dhp,dht,ir); // ! standard vapor
    rv0 = r; 
    ev0 = h - 1e5/r;
    wrsteamtab(1e5,300.0,0,r,drp,drt,h,dhp,dht,ir); // ! standard water
    rl0 = r; 
    el0 = h - 1e5/r;
    rg0 = 1.0;
    v0 = 1e-2*hc.vol;
    f10 = max(fabs(hc.etot), v0*rl0*el0);
    f20 = max(fabs(hc.rhov*hc.volg), v0*rv0);
    f30 = max(fabs(hc.rhol*hc.voll), v0*rl0);
    f40 = max(fabs(hc.rhon*hc.volg), v0*rg0);
    f50 = max(fabs(hc.rhoni*hc.volg), v0*rg0);
    f60 = hc.vol;

    // ! subcooled liq condition is externally specified (in massenevol)
    // !  p = psat(hc.tlc)
    // !  steamv(hc.tlc,p,v,h,s,cp,cv,hrho,cc,beta,akap,gam,0,i,0) // ! subcool water
    // !  hc.rholc = 1.0/v
    // !  hc.elc = h - p*v
    p = psatwaterf(hc.tlc);

    wrsteamtab(p,hc.tlc,0,r,drp,drt,h,dhp,dht,ir); // ! subcooled water
    hc.rholc = r; 
    hc.elc = h - p/r;
    hc.vollc = hc.mlc/hc.rholc;

    // ! initial guess of the evaporation/coondensation mass
    dmev = 0.0 ;

    // ! newton solution loop
    for (inew=1;inew<=nnewloop;inew++) // do inew = 1, nnewloop
    {
        // ! update thermophycical properties and related factors
        // ! vapor and gases
        p = hc.pv;
        // !     dtsp = dtsatdp(p)
        // !     steamv(ts,p,v,h,s,cp,cv,hrho,cc,beta,akap,gam,2,i,1) // ! sat steam
        // !     hc.rhov = 1.0/v; hc.ev = h - p*v; hc.cvv = cv; hc.cpv = cp
        ts = tsatwaterf(p);
        dtsp = dtsdpwaterf(p);
        // printf("getequil-11:%f, %d\n", hc.pv, inew);
        wrsteamtab(p,ts,1,r,drp,drt,h,dhp,dht,ir); // ! sat. steam
        hc.rhov = r; 
        hc.ev = h - p/r; 
        hc.cpv = dht; 
        hc.cvv = dht-ts/pow(r,2)*pow(drt,2)/drp;
        hc.ts = ts;
        hc.ccv = dht*hc.mmv/rgasgen; // ! ccv is calculated by precise cp
        // !     aevpv = (akap*p-beta*ts)*v + (cp-beta*p*v)*dtsp
        // !     arvpv = akap/v - beta/v*dtsp
        aevpv = dhp - 1.0/r + drp/pow(r,2)*p + (dht + drt/pow(r,2)*p)*dtsp;
        arvpv = drp + drt*dtsp;
        // printf("getequil-22:%f\n", hc.pv);

        stigas2(hc.pn,ts,hc.mmn,hc.ccn,hc.rhon,hc.en,hc.cvn,hc.cpn,drp,drt);
        aenpv = hc.cvn*dtsp;
        arnpv = drt*dtsp;
        arnpn = drp;
        stigas2(hc.pni,ts,hc.mmni,hc.ccni,hc.rhoni,hc.eni,hc.cvni,hc.cpni,drp,drt);
        aenipv = hc.cvni*dtsp;
        arnipv = drt*dtsp;
        arnipni = drp;

        hc.ptot = hc.pv + hc.pn + hc.pni;
        hc.rhotot = hc.rhov + hc.rhon + hc.rhoni;

        // ! water
        // ! consideration of condensation retard by non-condensible gases
        // ! if fnr=1 no effect
        hc.fnr = max(hc.fnr, 1.0);
        if((hc.fnr-1.0)*(hc.pn+hc.pni) > p-pliml) {hc.fnr = 1.0+(p-pliml)/(hc.pn+hc.pni);}
        p = p - (hc.fnr-1.0)*(hc.pn + hc.pni);
        // !     dtsp = dtsatdp(p)
        dtsp = dtsdpwaterf(p);
        if(itl == 0) { // ! tl is given by the equilibrium state
            // !  steamv(ts,p,v,h,s,cp,cv,hrho,cc,beta,akap,gam,4,i,1) // ! sat water
            // !  hc.tl = ts;
            // !  hc.rhol = 1.0/v; hc.el = h - p*v; hc.cpl = cp; hc.betal=beta
            // !  aelpv = (akap*p-beta*ts)*v + (cp-beta*p*v)*dtsp
            // !  arlpv = akap/v - beta/v*dtsp
            ts = tsatwaterf(p);
            wrsteamtab(p,ts,0, r,drp,drt,h,dhp,dht,ir); // ! sat. water
            hc.tl = ts;
            hc.rhol = r; 
            hc.el = h - p/r; 
            hc.cpl = dht; 
            hc.betal = -drt/r;
            aelpv = (drp*p+drt*ts)/pow(r,2) + (dht+drt/pow(r,2)*p)*dtsp;
            aelpn = -(hc.fnr-1.0)*aelpv;
            arlpv = drp + drt*dtsp;
            arlpn = -(hc.fnr-1.0)*arlpv;
        }
        else { // ! use externally given tl as a constant
            tl = hc.tl;
            // !  steamv(tl,p,v,h,s,cp,cv,hrho,cc,beta,akap,gam,0,i,1) // ! keep subcool
            // !  hc.rhol = 1.0/v; hc.el = h - p*v; hc.cpl = cp; hc.betal = beta
            // !  aelpv = (akap*p-beta*tl)*v
            // !  arlpv = akap/v 
            wrsteamtab(p,tl,0, r,drp,drt,h,dhp,dht,ir); // ! keep subcool
            hc.rhol = r; 
            hc.el = h - p/r; 
            hc.cpl = dht; 
            hc.betal = -drt/r;
            aelpv = (drp*p+drt*tl)/pow(r,2);
            arlpv = drp;
            aelpn = -(hc.fnr-1.0)*aelpv;
            arlpn = -(hc.fnr-1.0)*arlpv;
        }

        // ! present residuals
        f1 = (hc.mva + dmev)*hc.ev + hc.mn*hc.en + hc.mni*hc.eni + 
             (hc.mla - dmev)*hc.el + hc.mlc*hc.elc - hc.etot;
        f2 = hc.mva + dmev - hc.rhov*hc.volg;
        f3 = hc.mla - dmev - hc.rhol*hc.voll;
        f4 = hc.mn - hc.rhon*hc.volg;
        f5 = hc.mni - hc.rhoni*hc.volg;
        f6 = hc.volg + hc.voll + hc.vollc - hc.vol;

        res = fabs(f1/f10);
        res = max(res, fabs(f2/f20));
        res = max(res, fabs(f3/f30));
        res = max(res, fabs(f4/f40));
        res = max(res, fabs(f5/f50));
        res = max(res, fabs(f6/f60));

        if(res < epsnew) {break;}

        // ! jacobian matrix (factors on each corrections)
        // a = 0.0; // ! clear the jacobian area
        fill(a.begin(),a.end(),vector<long double>(7,0.0));

        // !         *primary var corrections
        // !*res_eq   ddmev      dpv        dpn        dpni       dvg      dvl 
        // ! df1      ev-el      *A         *B          *B
        // ! df2       1      -vg*arvpv                          -rhov
        // ! df3      -1      -vl*arlpv   -vl*arlpn  -vl*arlpn            -rhol
        // ! df4              -vg*arnpv   -vg*arnpn              -rhon
        // ! df5              -vg*arnipv             -vg*arnipni -rhoni
        // ! df6                                                   1        1
        // !       *A = (mva+dmev)*aevpv+mn*aenpv+mni*aenipv+(mla-dmev)*aelpv
        // !       *B = (mla-dmev)*aelpn
        a[0][0] = hc.ev-hc.el;
        a[0][1] = (hc.mva+dmev)*aevpv + hc.mn*aenpv + hc.mni*aenipv + (hc.mla-dmev)*aelpv;
        a[0][2] = (hc.mla-dmev)*aelpn; 
        a[0][3] = a[0][2];
        a[1][0] = 1.0; 
        a[1][1] = -hc.volg*arvpv; 
        a[1][4] = -hc.rhov;
        a[2][0] = -1.0; 
        a[2][1] = -hc.voll*arlpv; 
        a[2][2] = -hc.voll*arlpn; 
        a[2][3] = a[2][2];
        a[2][5] = -hc.rhol;
        a[3][1] = -hc.volg*arnpv; 
        a[3][2] = -hc.volg*arnpn; 
        a[3][4] = -hc.rhon;
        a[4][1] = -hc.volg*arnipv; 
        a[4][3] = -hc.volg*arnipni; 
        a[4][4] = -hc.rhoni;
        a[5][4] = 1.0; 
        a[5][5] = 1.0;

        // ! res_eqs to be solved : dfi = -fi][ i=1--6
        a[0][6] = -f1;
        a[1][6] = -f2;
        a[2][6] = -f3;
        a[3][6] = -f4;
        a[4][6] = -f5;
        a[5][6] = -f6;

        gauss(a, 6, is, neps, 1); // ! is = 1 if fail in gauss
        ddmev = a[0][6]; 
        dpv = a[1][6]; 
        dpn = a[2][6]; 
        dpni = a[3][6]; 
        dvg = a[4][6]; 
        dvl = a[5][6];


        // printf("getequil-22:%f, clear\n", hc.pv);
        // printf("again %f\n", dpv);

        dmev = dmev + ddmev;
        hc.dmwev = dmev;
        hc.pv = hc.pv + dpv ; 
        hc.pv = max(hc.pv, pliml);
        hc.pn = hc.pn + dpn ; 
        hc.pn = max(hc.pn, 0.0);
        hc.pni = hc.pni + dpni ; 
        hc.pni = max(hc.pni, 0.0);
        hc.volg = hc.volg + dvg; 
        hc.volg = max(hc.volg, 0.0);
        hc.voll = hc.voll + dvl; 
        hc.voll = max(hc.voll, 0.0);
        hc.mv = hc.mva + dmev ; 
        hc.mv = max(hc.mv, 0.0);
        hc.ml = hc.mla - dmev ; 
        hc.ml = max(hc.ml, 0.0);

        // printf("getequil-33:%f\n", hc.pv);

    }

    // ! if ml <= 0 and dmev >0, solve equations for superheated vapor
    if(res > epsnew && hc.ml <=0.0 && dmev > 0.0) {
        hc.pv = pvorg; 
        hc.pn = pnorg; 
        hc.pni = pniorg; 
        hc.ts = tsorg;

        // ! no water, no phase change
        hc.volg = hc.vol - hc.vollc;
        hc.voll = 0.0; 
        hc.ml = 0.0;
        hc.dmwev = 0.0;
        hc.mv = max(hc.mva + hc.mla, 0.0); // ! in-flow liquid should be added
                                        // ! vaporized because ml<0 in the normal scheme
        if(hc.volg <= 0.0) {
            // errorhdr("getequil: volg <= 0// !", "stop");
            printf("getequil: volg <= 0// !\n");
            exit(EXIT_FAILURE);
        }

        // !  hc.rhov = hc.mv/hc.volg
        // !  hc.rhon = hc.mn/hc.volg
        // !  hc.rhoni = hc.mni/hc.volg

        // ! initial guess of the vapor state = (hc.pv, hc.ts)
        // ! consider the stability limit of the vapor state
        ts1 = tsatwaterf(hc.pv);
        if (hc.ts < ts1) {hc.ts = ts1;}

        // ! newton loop for super heated vapor and gas
        for(inew=1; inew<=nnewloop;inew++) // do inew = 1, nnewloop
        {
            p = hc.pv; 
            ts = hc.ts;
            wrsteamtab(p,ts,1,r,drp,drt,h,dhp,dht,ir); // ! steam, can be superheated
            hc.rhov = r; 
            hc.ev = h - p/r; 
            hc.cpv = dht; 
            hc.cvv = dht-ts/pow(r,2)*pow(drt,2)/drp;
            hc.ccv = dht*hc.mmv/rgasgen; // ! ccv is calculated by precise cp
            aevpv = dhp - 1.0/r + drp/pow(r,2)*p; 
            aevt =  dht + drt/pow(r,2)*p;
            arvpv = drp; 
            arvt = drt;

            stigas2(hc.pn,ts,hc.mmn,hc.ccn,hc.rhon,hc.en,hc.cvn,hc.cpn,drp,drt);
            aent = hc.cvn; // ! for ideal gas (dedt)(p=const)=cv
            arnpn = drp; arnt = drt;
            stigas2(hc.pni,ts,hc.mmni,hc.ccni,hc.rhoni,hc.eni,hc.cvni,hc.cpni,drp,drt);
            aenit = hc.cvni;
            arnipni = drp; arnit = drt;

            hc.ptot = hc.pv + hc.pn + hc.pni;
            hc.rhotot = hc.rhov + hc.rhon + hc.rhoni;

            // ! present residuals
            f1 = hc.mv*hc.ev + hc.mn*hc.en + hc.mni*hc.eni + hc.mlc*hc.elc - hc.etot;
            f2 = hc.mv - hc.rhov*hc.volg;
            f4 = hc.mn - hc.rhon*hc.volg;
            f5 = hc.mni - hc.rhoni*hc.volg;

            res = fabs(f1/f10);
            res = max(res, fabs(f2/f20));
            res = max(res, fabs(f4/f40));
            res = max(res, fabs(f5/f50));
    
            if(res < epsnew) {break;}

            // ! jacobian matrix (factors on each corrections)
            // a1 = 0.0 // ! clear the jacobian area
            fill(a1.begin(),a1.end(),vector<long double>(5,0.0));

            // !         *primary var corrections
            // !*res_eq   dts         dpv        dpn        dpni  
            // ! df1      *A       mv*aevpv
            // ! df2    -vg*arvt  -vg*arvpv
            // ! df4    -vg*arnt              -vg*arnpn
            // ! df5    -vg*arnit                       -vg*arnipni
            // !     *A = mv*aevt+mn*aent+mni*aenit
            a1[0][0] = hc.mv*aevt + hc.mn*aent + hc.mni*aenit;
            a1[0][1] = hc.mv*aevpv;
            a1[1][0] = -hc.volg*arvt; 
            a1[1][1] = -hc.volg*arvpv;
            a1[2][0] = -hc.volg*arnt; 
            a1[2][2] = -hc.volg*arnpn;
            a1[3][0] = -hc.volg*arnit; 
            a1[3][3] = -hc.volg*arnipni;

            // ! res_eqs to be solved : dfi = -fi][ i=1--4
            a1[0][4] = -f1;
            a1[1][4] = -f2;
            a1[2][4] = -f4;
            a1[3][4] = -f5;

            gauss(a1, 4, is, neps, 1); // ! is = 1 if fail in gauss
            dts = a1[0][4]; 
            dpv = a1[1][4]; 
            dpn = a1[2][4]; 
            dpni = a1[3][4];

            hc.ts = hc.ts + dts ; 
            hc.ts = max(hc.ts, tliml);
            hc.pv = hc.pv + dpv ; 
            hc.pv = max(hc.pv, pliml);
            hc.pn = hc.pn + dpn ; 
            hc.pn = max(hc.pn, 0.0);
            hc.pni = hc.pni + dpni ; 
            hc.pni = max(hc.pni, 0.0);
        }
        f3 = 0.0; 
        f6 = 0.0;
    }

    // ! show residuals and stop if not converged
    if(res > epsnew) {
        // write(*,"(a)") "# getequil: Newton loop NOT converged. Residuals are:"
        // write(*,"(a,1pe12.4,a,1pe12.4)") "#  f1:",f1," f1/f10: ",f1/f10
        // write(*,"(a,1pe12.4,a,1pe12.4)") "#  f2:",f2," f2/f20: ",f2/f20
        // write(*,"(a,1pe12.4,a,1pe12.4)") "#  f3:",f3," f3/f30: ",f3/f30
        // write(*,"(a,1pe12.4,a,1pe12.4)") "#  f4:",f4," f4/f40: ",f4/f40
        // write(*,"(a,1pe12.4,a,1pe12.4)") "#  f5:",f5," f5/f50: ",f5/f50
        // write(*,"(a,1pe12.4,a,1pe12.4)") "#  f6:",f6," f6/f60: ",f6/f60
        // write(*,"(a,1pe12.4)") "#  Allowed limit: epsnew=", epsnew
        printf("# getequil: Newton loop NOT converged. Residuals are:\n");
        printf("#  f1: %f, f1/f10: %f \n", f1, f1/f10);
        printf("#  f2: %f, f2/f20: %f \n", f2, f2/f20);
        printf("#  f3: %f, f3/f30: %f \n", f3, f3/f30);
        printf("#  f4: %f, f4/f40: %f \n", f4, f4/f40);
        printf("#  f5: %f, f5/f50: %f \n", f5, f5/f50);
        printf("#  f6: %f, f6/f60: %f \n", f6, f6/f60);
        printf("#  Allowed limit: epsnew= %f\n", epsnew);
        is = -1; // ! is = -1 if fail in newton
    }

    // ! check overshooting of pressure/temperature for dt control
    if(hc.mdov+hc.mdon+hc.mdoni+hc.mdol+hc.mdorcic+hc.mdolc>0.0 &&
        hc.ptot < ptotorg && hc.ptot < hc.pout*0.99 && is == 0) {
        is = -2; // ! is=-2 if overshooting pout with successful solution
        }
    if(hc.qsrcex != 0.0 && (tsorg - hc.tsrcex)*(hc.ts - hc.tsrcex) < 0.0 
        && fabs(hc.ts - hc.tsrcex) > 2.0) {is = -3;}
    // ! do the water level
    if(hc.tvollev.n>1) {
        itabeval(hc.tvollev, hc.vollc, hc.hlc) ;
        itabeval(hc.tvollev, hc.voll+hc.vollc, hc.hl) ;
        itabeval(hc.tvollev, hc.vol, hc.htop) ;
    }

    ist = is;

}
/*###############################################################################*/

// int main ()
// {
//     double test10[10]={1.0,};
//     vector<double> test20(10,1.0);
//     htcan hc1;

//     cout<< test20[1] << endl;
//     printf("TEST\n");
//     printf("FIN\n");
//     return 0;
// }

#endif