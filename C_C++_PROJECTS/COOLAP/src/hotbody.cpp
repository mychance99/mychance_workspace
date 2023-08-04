#ifndef HOTBODY_H
#define HOTBODY_H

#include <iostream>
#include <vector>
#include <string>
// #include <cmath>
#include <math.h>

using namespace std;


#include "inttable_mod.cpp"
#include "coolpropdef.cpp"
#include "heattr.cpp"
#include "consts.cpp"
#include "idgas.cpp"
#include "math_mod.cpp"

/*###############################################################################*/ 

const int nmaxqgen=10;
double temp0=273.15;

struct htbody
{
    int ihtr=1,   // ! switch for the heat transfer model 
            // ! 1:use (calc. by correlations), 0:not use (give htc as constants 
            // !   htcwet, htcdry)
        icdn=0,   // ! switch for condensation heat transfer in dry zone
            // ! 0:not use (default), 1:use (nusselt laminar film condensation correlation
            // !   with modification factor fhtcdn)
        irsf=1,   // ! switch for the heat transfer efficiency of the dry surface, rsfdr
            // ! 1:use (HTC at dry zone attenuated), 0:not use
        iasf=1,   // ! switch for the surface area handling
            // ! 1:give asf or tasf, 0:give ll(=d/2) and hh and evaluate asf by
            // ! asf=vol*(4/d+2/h) according to bundle model
        nqgen=0;  // ! number of heat generation sources (e.g. different mechanisms, etc.)

    double  hlow=0.0, hhigh=0.0,  // ! lower/upper end position of the body (m)
                       // ! corresponding with "tvollev" of the can.
            hh=0.0, ll=0.0,   // ! height (hhigh-hlow) (m) and
                       // ! thermal boundary layer size inside the body 
                       // ! e.g. radius, thickness/2
            acr=0.0, asf=0.0, // ! cross section (horizontal) area and surface area (m2)
            vol,rsfdr,  // ! volume of the body (m3), dry zone heat transfer surface efficiency,
            fwet=0.0,     // ! wet zone height (under water) (m)
            twet, tdry, twsf, tdsf,  // ! temperatures(K) (wet, dry zones, and their surface)
            tini=0.0,  // ! initial temperature for input
            rho, cp, lam, epsrad=0.0,  // ! properties of the body 
                      // ! (density(kg/m3), specific heat(J/kgK), heat conductivity(W/mK), emmissivity)
            qgen=0.0,qtr=0.0,qrlx=0.0, // ! powers (W) (generation, transfer to the can, between the dry  wet zones)
            htcwet=0.0,htcdry=0.0,  // ! given constant heat transfer coefficients (W/m2K)
            qflxwet, qflxdry,  // ! surface heat fluxes [W/m2] in the wet and dry zones
                // ! direction of heat flow positive: body=>coolant, dry zone=>wet zone
            ftdrp=0.5,  // ! tav-tsf=(q"*ll/2/lam)*ftdrp 
                           // ! 1/4 for cylinders, 2/3 for flat plates
            fhtr=1.0,  // ! factor for heat transfer to the coolant
            frlx=2.0,  // ! factor for internal heat flow dry zone=>wet zone
            ene=0.0,intqgen=0.0, // ! energy (J) and time integral of the generated heat (J)
            fhtdrywet=0.0, // ! switch for the usage of wet zone HTC for the dry zone
                        // ! >0:use the wet HTC for the dry zone with this factor, 0:no (use the dry zone HTC)
            fhtcdn=1.0;  // ! factor for condensation heat transfer correlation for dry zone
        vector<double> ftqgen = vector<double>(nmaxqgen,1.0);  // ! a factor to modify tqgen()
    string  htmodw,  // ! heat transfer mode of the wet zone
                    // ! ([c]onv./[n]ucl. boil./[f]ilm boil./[t]ransition boil./- no wet zone)
            htmoddw; // ! heat transfer mode of the dry zone where wet zone HTC is applied
    timetab tqgen[nmaxqgen],  // ! time table for qgen (if given as time dependent
            tftqgen[nmaxqgen], // ! constants), and modification factor for tqgen => set ftqgen
            tfhtr,  // !  time tables for fhtr
            tfhtdrywet,  // !  time tables for fhtdrywet
            thlow, thhigh, tasf,tll, // ! time tables for hlow,hhigh,asf,ll
                  // !  => movement and deformation
            tacr; // ! time table for acr; if not given, acr is calculated by acr=vol/(hhigh-hlow)
                  // ! so that the initial volume is preserved; if given, the volume (mass) can change
    coolprop clp;
};
/*###############################################################################*/ 







/*###############################################################################*/ 
void initbody(htbody &hb)
{
    int ihtr, i;
    double hlow, hhigh, asf, ll,acr,vol,tini,rho,cp,lam,epsrad;

    for (i=1; i<=hb.nqgen; i++)
    {
        if(hb.tqgen[i].n > 0)  {ttabcheck(hb.tqgen[i]);}
        if(hb.tftqgen[i].n > 0)  {ttabcheck(hb.tftqgen[i]);}
    }

    if(hb.tfhtr.n > 0)  {ttabcheck(hb.tfhtr);}
    if(hb.thlow.n > 0)  {ttabcheck(hb.thlow);}
    if(hb.thhigh.n > 0)  {ttabcheck(hb.thhigh);}
    if(hb.tasf.n > 0)  {ttabcheck(hb.tasf);}
    if(hb.tacr.n > 0)  {ttabcheck(hb.tacr);}
    if(hb.tll.n > 0)  {ttabcheck(hb.tll);}
    if(hb.tfhtdrywet.n > 0)  {ttabcheck(hb.tfhtdrywet);}

    // ! geometry etc.
    hlow = hb.hlow; hhigh = hb.hhigh;
    hb.hh = hb.hhigh - hb.hlow;
    if(hb.hh <= 0.0)  {
        printf("TMP-hotbody-initbody\n");
        exit(EXIT_FAILURE);
    }// errorhdr("initbody: hhigh < hlow// !", "cont")
    tini = hb.tini;
    hb.twet = tini; hb.tdry = tini;
    hb.twsf = tini; hb.tdsf = tini;
    // !  rho = hb.rho; cp = hb.cp; lam = hb.lam; epsrad = hb.epsrad
    if(hb.vol <= 0.0) {hb.vol = hb.acr*hb.hh;}
    if(hb.acr <= 0.0) {hb.acr = hb.vol/hb.hh;}
    // !  ll = hb.ll; acr = hb.acr; vol = hb.vol
    if(hb.vol <= 0.0 || hb.acr <= 0.0) {
        //  errorhdr("initbody: zero acr or vol// !","cont")
        printf("TMP-hotbody-initbody\n");
        exit(EXIT_FAILURE);
    }
    if(hb.tini <= 0.0) { // errorhdr("initbody: invalid tini// !","cont")
        printf("TMP-hotbody-initbody\n");
        exit(EXIT_FAILURE);
    }

    hb.htmodw = '-'; hb.htmoddw = '-';

    // ! bundle geometry: pin dia.=d=2*ll, number=n
    // !   asf=pi*d*h*n, acr=pi/4*d^2*n, vol=acr*h
    // !   asf/vol=4/d+2/h=2/ll+2/h
    if(hb.iasf != 1) {
        hb.asf = hb.vol*(2.0/hb.ll + 2.0/hb.hh);
    }
}
/*###############################################################################*/ 

/*###############################################################################*/
void htcwet(double time,double dt,double kb,htbody &hb,int &is,string &cmode,double &hwet)
// ! calculation of the heat transfer coefficient
{
    double tspt,hfg,sig,tav,tc,tsf,tm,betavf,rhovf,
    drp,drt,h,dhp,dht,cpvf,muvf,lamvf,rhov,qflux,qboil,qcon,qrad,
    tbi,tchf,tmfb,qchf,qmfb,q0,rq0,t0,rq1,t1,kbt,qbc,qbm,rqm;
    int i,ir;
    cmode = "-"; 
    is = 0;

    // ! collect phys. properties
    hb.clp.mul = viscwaterf(hb.clp.rhol, hb.clp.tl);
    hb.clp.laml = tconwaterf(hb.clp.rhol, hb.clp.tl);
    // !  call steamz(tspt,hb.clp.ptot,v,hfg,s,2) // ! sat.vapor p base
    // !  rhov = 1.0/v
    tspt = tsatwaterf(hb.clp.ptot);
    wrsteamtab(hb.clp.ptot,tspt,1,rhov,drp,drt,h,dhp,dht,ir); // ! sat.vapor
    // !  call steamz(tspt,hb.clp.ptot,v,h,s,4) // ! sat.water p base
    // !  hfg = hfg - h
    hfg = hfgwaterf(hb.clp.ptot);
    sig = stenwaterf(tspt);

    tav = hb.twet;
    tc = hb.clp.tl;
    tsf = hb.twsf; // ! the previous value as initial guess
    if(fabs(tav-tsf) > fabs(tav-tc)) {tsf = 0.5*(tav+tc);} // ! bounding

    // ! vapor film properties 
    tm = (tspt + max(tspt+1.0, tsf))/2.0;
    betavf = 1.0/tm;
    // !  call steam(tm,hb.clp.ptot,v,h,s,cpvf,0) // ! super heated vapor
    // !  rhovf = 1.0/v
    wrsteamtab(hb.clp.ptot,tm,1,rhovf,drp,drt,h,dhp,cpvf,ir); // ! super heated vapor
    muvf = viscwaterf(rhovf, tm);
    lamvf = tconwaterf(rhovf, tm);

    t0 = tsf;
    htransfl(t0, hb.epsrad, hb.clp.ptot, tspt, hfg,
        tc,hb.clp.rhol,hb.clp.cpl,hb.clp.mul,hb.clp.laml,sig,hb.clp.betal,
        rhov,rhovf,cpvf,muvf,lamvf,betavf, 
        qflux, qboil, qcon, qrad, cmode, tbi,tchf,tmfb,qchf,qmfb);
    q0 = hb.fhtr*qflux;
    rq0 = q0 + kb*(t0 - tav);  // ! the function to be zeroed

    // ! -- bi-section method to get tsf ---

    // ! set the bounding temperatures: t0 and t1
    if(cmode != "c") { // ! boiling
        kbt = -(qchf-qmfb)/(tchf-tmfb); // ! gradient of trans. boil. curve
        qbc = -kb*(tchf-tav);   // ! qb at tchf
        qbm = -kb*(tmfb-tav); // ! qb at tminfb
        if (kbt < kb) {
            // ! stable single solution in trans. zone is possible
            // ! all zones are uniform in view of getting the solution by bi-section
            if( rq0 > 0.0) { // ! solution is in the LEFT
                t1 = min(tc, tav);
            }
            else{ // ! solution is in the RIGHT
                t1 = max(tav, tc);
            }
        }
        else{
            // ! transition zone is unstable
            // ! leap "n" <-> "f" should be considered
            if (cmode == "t") { // ! do the transition first.
                if(qbc <= qchf) { // ! leap to nucl.
                    t0 = tchf; 
                    t1 = min(tc, tav);
                    rq0 = 0.0;
                }
                else { // ! leap to film boil.
                    t0 = tmfb; 
                    t1 = max(tav, tc);
                    rq0 = 0.0;
                }
            }
            else {
                if( rq0 > 0.0) { // ! solution is in the LEFT
                    if(cmode == "n") {
                        t1 = min(tc, tav);
                    }
                    else { // ! "f"
                        if(qbm >= qmfb) { // ! within film boiling
                            t1 = tmfb;
                        }
                        else { // ! leap to nucl. boil or conv.
                            t0 = tchf; 
                            t1 = min(tc, tav);
                            rq0 = 0.0;
                        }
                    }
                }
                else { // ! solution is in the RIGHT
                    if (cmode == "n") {
                        if(qbc <= qchf) { // ! nucl.
                            t1 = tchf;
                        }
                        else { // ! leap to film boil.
                            t0 = tmfb; 
                            t1 = max(tav, tc);
                            rq0 = 0.0;
                        }
                    }
                    else { // ! "f"
                        t1 = max(tav, tc);
                    }
                }
            }
        }
    }
    else { // ! convection
        if( rq0 > 0.0) { // ! solution is in the LEFT
            t1 = min(tc, tav);
            for (i=1; i<=50; i++) // do i=1,50
            {
                htransfl(t1, hb.epsrad, hb.clp.ptot, tspt, hfg,
                    tc,hb.clp.rhol,hb.clp.cpl,hb.clp.mul,hb.clp.laml,sig,hb.clp.betal,
                    rhov,rhovf,cpvf,muvf,lamvf,betavf, 
                    qflux, qboil, qcon, qrad, cmode, tbi,tchf,tmfb,qchf,qmfb);
                if (hb.fhtr*qflux < -kb*(t1 - tav)) {break;}// exit
                t1 = t1 - 10.0;
            }
        }
        else {
            t1 = t0 + 10.0;
            for (i=1; i<=50; i++) // do i=1,50
            {
                // !           t1 = min(t1, max(tav,tc))
                htransfl(t1, hb.epsrad, hb.clp.ptot, tspt, hfg,
                    tc,hb.clp.rhol,hb.clp.cpl,hb.clp.mul,hb.clp.laml,sig,hb.clp.betal,
                    rhov,rhovf,cpvf,muvf,lamvf,betavf, 
                    qflux, qboil, qcon, qrad, cmode, tbi,tchf,tmfb,qchf,qmfb);
                if (hb.fhtr*qflux > -kb*(t1 - tav)) {break;} // exit
                t1 = t1 + 10.0;
            }
        }
    }

    if(rq0 == 0.0) {
        htransfl(t0, hb.epsrad, hb.clp.ptot, tspt, hfg,
        tc,hb.clp.rhol,hb.clp.cpl,hb.clp.mul,hb.clp.laml,sig,hb.clp.betal,
        rhov,rhovf,cpvf,muvf,lamvf,betavf, 
        qflux, qboil, qcon, qrad, cmode, tbi,tchf,tmfb,qchf,qmfb);
        rq0 = hb.fhtr*qflux + kb*(t0 - tav);
    }
    htransfl(t1, hb.epsrad, hb.clp.ptot, tspt, hfg,
        tc,hb.clp.rhol,hb.clp.cpl,hb.clp.mul,hb.clp.laml,sig,hb.clp.betal,
        rhov,rhovf,cpvf,muvf,lamvf,betavf, 
        qflux, qboil, qcon, qrad, cmode, tbi,tchf,tmfb,qchf,qmfb);
    rq1 = hb.fhtr*qflux + kb*(t1 - tav);

    // ! bi-section convergence
    is = 1;
    for(i=1; i<=20; i++) // do i = 1,20
    {
        tm = 0.5*(t0+t1);
        htransfl(tm,hb.epsrad, hb.clp.ptot, tspt, hfg,
            tc,hb.clp.rhol,hb.clp.cpl,hb.clp.mul,hb.clp.laml,sig,hb.clp.betal,
            rhov,rhovf,cpvf,muvf,lamvf,betavf, 
            qflux, qboil, qcon, qrad, cmode, tbi,tchf,tmfb,qchf,qmfb);
        rqm = hb.fhtr*qflux + kb*(tm - tav);
        if(rq0*rqm > 0.0) { // ! solution is in m--1
            t0 = tm; 
            rq0 = rqm;
        }
        else {             // ! solution is in 0--m
            t1 = tm; 
            rq1 = rqm;
        }
        if(fabs(t0 - t1) < 1e-2) {
            is = 0; 
            break; // exit
        }
        if(rq1*rq0 > 0.0) {
            // write(*,'(a,1p,2(1x,e11.3))') &
            // 'hbevol: bi-section for tsf (wet) failed. time,dt=',time,dt
            printf("hbevol: bi-section for tsf (wet) failed. time,dt=%f %f", time, dt);
            // write(*,'(a,1p,7(1x,e10.2))') &
            // '   t0,t1,tav,tc,qflux,rq0,rq1=',t0,t1,tav,tc,qflux,rq0,rq1
            printf("   t0,t1,tav,tc,qflux,rq0,rq1=%f %f %f %f %f %f %f",t0,t1,tav,tc,qflux,rq0,rq1);
            break; // exit
        }
    }
    tsf = 0.5*(t0 + t1);
    // ! -- bi-section method to get tsf --end--
    hb.htmodw = cmode;
    if(fabs(tsf - tc) > 0.0) {
        hwet = fabs(hb.fhtr*qflux/(tsf - tc));
    }
    else {hwet = 0.0;}
    // ! only the heat transfer coeff. is obtained here.
    // ! qflux = htc * (tsf - tc)
}
/*###############################################################################*/




/*###############################################################################*/
void htcdry (double time,double dt,double kb,htbody &hb,int &is,double &hdry)
{
    double tav,tc,tsf,tm,betagf,rhogf,
    cpgf,mugf,lamgf,qflux,qboil,qcon,qrad,
    tbi,tchf,tmfb,qchf,qmfb,q0,rq0,t0,rq1,t1,kbt,qbc,qbm,rqm;
    int i,ir;

    is = 0;

    printf("#############################################\n");
    printf("     not coded...    hotbody.htcdry          \n");
    printf("#############################################\n");
/* 
    // ! collect phys. properties (no bulk property is referred below at present)
    mixpropgas(hb.clp.tg,hb.clp.ng,hb.clp.gname,hb.clp.pg, 
        hb.clp.mug,hb.clp.lamg,hb.clp.cpg,hb.clp.betag);

    tav = hb.tdry;
    tc = hb.clp.tg;
    tsf = hb.tdsf; // ! the previous value as initial guess
    if(abs(tav-tsf) > abs(tav-tc)) {tsf = 0.5*(tav+tc);} // ! bounding

    // ! heat transfer surface efficiency in the dry zone
    // ! (assume a big rod of cross section acr (instead of a bundle with more surface area),
    // !  height hh; especially with radiation, the surface does not act effectively)
    if(hb.irsf == 1) {
        if(hb.asf <= 0.0) // call errorhdr("hbevol: zero or less surface area// !","stop")
        {
            printf("hbevol: zero or less surface area// !\n");
            exit(EXIT_FAILURE);
        }
        hb.rsfdr = (sqrt(4.0*pi*hb.acr)*hb.hh+2.0*hb.acr)/hb.asf;
    }
    else {
        hb.rsfdr = 1.0;
    }

    // ! gas film properties
    tm = 0.5*(tc + tsf);
    mixpropgas(tm,hb.clp.ng,hb.clp.gname,hb.clp.pg,
        mugf,lamgf,cpgf,betagf);
    rhogf = hb.clp.rhog * tc/tm;

    // ! bi-section to get tsf
    t0 = tsf;
    htransfg(t0, hb.epsrad, tc, rhogf,cpgf,mugf,lamgf,betagf, 
        qflux, qcon, qrad);
    rq0 = hb.fhtr*qflux + kb*(t0 - tav); // ! residual of q" = func. that should be 0

    if(rq0 > 0.0) {// ! solution is in the LEFT
        t1 = min(tc, tav);
    }
    else { // ! in the RIGHT
        t1 = max(tav, tc);
    }
    htransfg(t1, hb.epsrad, tc, rhogf,cpgf,mugf,lamgf,betagf, 
        qflux, qcon, qrad);
    rq1 = hb.fhtr*qflux + kb*(t1 - tav);

    // ! bi-section converge
    is = 1;
    for (i=1; i<=20; i++) // do i=1,20
    {
        tm = 0.5*(t0+t1);
        htransfg(tm, hb.epsrad, tc, rhogf,cpgf,mugf,lamgf,betagf, 
            qflux, qcon, qrad);
        rqm = hb.fhtr*qflux + kb*(tm - tav);
        if(rq0*rqm > 0.0) { // ! solution is in m--1
            t0 = tm; rq0 = rqm;
        }
        else {            // ! solution is in 0--m
            t1 = tm; rq1 = rqm;
        }
        if(abs(t0 - t1) < 1e-2) {
            is = 0; exit(EXIT_FAILURE);
        }
        if(rq1*rq0 > 0.0) {
            // write(*,'(a,1p,2(1x,e11.3))') &
            // 'hbevol: bi-section for tsf (dry) failed. time,dt=',time,dt
            // write(*,'(a,1p,7(1x,e10.2))') &
            // '   t0,t1,tav,tc,qflux,rq0,rq1=',t0,t1,tav,tc,qflux,rq0,rq1
            printf("hbevol: bi-section for tsf (dry) failed. time,dt= %f %f",time,dt);
            printf("   t0,t1,tav,tc,qflux,rq0,rq1=%f %f %f %f %f %f %f",t0,t1,tav,tc,qflux,rq0,rq1);
            exit(EXIT_FAILURE);
        }
    }
    tsf = 0.5*(t0 + t1);
    // ! -- bi-section method to get tsf --end--
    // !  hb.qflxdry = qflux; hb.tdsf = tsf
    if(abs(tsf - tc) > 0.0) {
        hdry = abs(hb.fhtr*qflux/(tsf - tc)) * hb.rsfdr;
    }
    else {
        hdry = 0.0;
    }
    // ! only the heat transfer coeff. is obtained here.
    // ! qflux = htc * (tsf - tc)
    // ! reduction of HTC by surface efficiency is included (rsfdr)
     */
}

/*###############################################################################*/


/*###############################################################################*/
void htcdn(double tsf,double tsat,double hh,double grav,double rhol,double laml,
           double mul,double dhfg, 
           double &hcdn,double &qcdn)
// ! condensation heat transfer (with natural convection on a vertical wall)
{
    // ! nusselt condensation model (laminar, as lower end condensation htc)
    // ! nu=h*L/laml=0.943*(g*rhol*(rhol-rhog)*dhfg*L**3/mul/laml/(tsat-tw))**0.25
    // ! h=0.943*(g*rhol*(rhol-rhog)*dhfg/mul/(tsat-tw))**0.25 * laml**0.75 / L**0.25
    // ! q=h*(tsat-tw)=0.943*(g*rhol*(rhol-rhog)*dhfg/mul)**0.25 * ((tsat-tw)*laml)**0.75 / L**0.25    
    qcdn=-0.943*pow((grav*pow(rhol,2)/mul*dhfg),0.25) * 
        pow((max(tsat-tsf,0.0)*laml),0.75) / pow(max(hh,0.1),0.25);
    hcdn=fabs(qcdn/max(tsat-tsf,1.0));
}
/*###############################################################################*/


/*###############################################################################*/
void htcdrycdn (double time,double dt,double kb,htbody &hb,int &is,double &hdry)
// ! dry zone condensation heat transfer
// ! nusselt laminar film condensation correlation
// ! presently, no consideration on non-condensible gas effect (heat resistance by 
// ! accumulation near the interface)
{
    double tav,tc,tsf,tm,rhol,laml,mul,dhfg,
        qflux,htc, 
        tbi,tchf,tmfb,qchf,qmfb,q0,rq0,t0,rq1,t1,kbt,qbc,qbm,rqm;
    int i,ir;

    is = 0;

    tav = hb.tdry;
    tc = hb.clp.tg; // ! tg=tc=tl=sat. at steam partial pressure
    tsf = hb.tdsf; // ! the previous value as initial guess
    if(fabs(tav-tsf) > fabs(tav-tc)) {tsf = 0.5*(tav+tc);} // ! bounding

    // ! discard if the condition is not adequate
    if(tsf > tc) {
        hdry = 0.0;
        return;
    }

    // ! collect liquid phys. props.
    hb.clp.mul = viscwaterf(hb.clp.rhol, hb.clp.tl);
    hb.clp.laml = tconwaterf(hb.clp.rhol, hb.clp.tl);
    dhfg = hfgwaterf(hb.clp.pg[0]);

    // ! bi-section to get tsf
    t0 = tsf;
    htcdn(t0,tc,hb.hh,grav, hb.clp.cpl,hb.clp.laml,hb.clp.mul,dhfg,
          htc,qflux);
    rq0 = hb.fhtcdn*qflux + kb*(t0 - tav); // ! residual of q" = func. that should be 0

    if(rq0 > 0.0) { // ! solution is in the LEFT
        t1 = min(tc, tav);
    }
    else { // ! in the RIGHT
        t1 = max(tav, tc);
    }
    htcdn(t1,tc,hb.hh,grav, hb.clp.cpl,hb.clp.laml,hb.clp.mul,dhfg,
        htc,qflux);
    rq1 = hb.fhtcdn*qflux + kb*(t1 - tav);

    // ! bi-section converge
    is = 1;
    for (i=1; i<=20; i++) // do i=1,20
    {
        tm = 0.5*(t0+t1);
        htcdn(tm,tc,hb.hh,grav, hb.clp.cpl,hb.clp.laml,hb.clp.mul,dhfg,
            htc,qflux);
        rqm = hb.fhtcdn*qflux + kb*(tm - tav);
        if(rq0*rqm > 0.0) {// ! solution is in m--1
            t0 = tm; 
            rq0 = rqm;
        }
        else {            // ! solution is in 0--m
            t1 = tm; 
            rq1 = rqm;
        }
        if(fabs(t0 - t1) < 1e-2) {
            is = 0; 
            break;
        }
        if(rq1*rq0 > 0.0) {
            // write(*,'(a,1p,2(1x,e11.3))') &
            // 'hbevol: bi-section for tsf (drycdn) failed. time,dt=',time,dt
            // write(*,'(a,1p,7(1x,e10.2))') &
            // '   t0,t1,tav,tc,qflux,rq0,rq1=',t0,t1,tav,tc,qflux,rq0,rq1
            printf("hbevol: bi-section for tsf (drycdn) failed. time,dt= %f %f ",time,dt);
            printf("   t0,t1,tav,tc,qflux,rq0,rq1=%f %f %f %f %f %f %f",t0,t1,tav,tc,qflux,rq0,rq1);
            break;
        }
    }
    tsf = 0.5*(t0 + t1);
    // ! -- bi-section method to get tsf --end--

    hdry = hb.fhtcdn*htc;
}
/*###############################################################################*/



























/*###############################################################################*/ 
void hbevol (htbody &hb, double time, double dt, int ist=1)
// ! calc. heat transfer (set qsrc)
// ! calc. energy evolution of a body
// ! requitement before calling this
// !   set hb%fwet
// !   set coolant states/properties in hb%clp
// !     (ptot,ng,gname(),pg(),tg,tl,rhog(mixture),rhol,cpl,betal)
{
    int i, is1, is2, is3, neps, ir;
    double lrlx, x1, x2, x3, hwet,hdry,hdrywet,q0,q1,kb;
    // double a[3][2];
    string cmode;
    double dtemplim = 5.0;
    vector<vector<long double>> a = 
        vector<vector<long double>>(2, vector<long double>(3,0));


    // if(present(ist)) ist = 0
    if (ist == 1){ist = 0;}

    // ! do the time tabled conditions
    // ! heat soruce and transfer
    q0 = 0.0;
    for(i=0; i<hb.nqgen; i++) // do i = 1, hb.nqgen
    {
        q1 = 0.0;
        if(hb.tqgen[i].n > 0) {ttabeval(hb.tqgen[i], time, q1);}
        if(hb.tftqgen[i].n > 0) {ttabeval(hb.tftqgen[i], time, hb.ftqgen[i]);}
        q0 = q0 + q1 * hb.ftqgen[i];
    }
    // end do
    hb.qgen = q0;
    if(hb.tfhtr.n > 0) {ttabeval(hb.tfhtr, time, hb.fhtr);}
    // ! movement and deformation
    if(hb.thlow.n > 0) {ttabeval(hb.thlow, time, hb.hlow);}
    if(hb.thhigh.n > 0) {ttabeval(hb.thhigh, time, hb.hhigh);}
    if(hb.tasf.n > 0) {ttabeval(hb.tasf, time, hb.asf);}
    if(hb.tacr.n > 0) {ttabeval(hb.tacr, time, hb.acr);}
    if(hb.tll.n > 0) {ttabeval(hb.tll, time, hb.ll);}
    if(hb.tfhtdrywet.n > 0) {ttabeval(hb.tfhtdrywet, time, hb.fhtdrywet);}
    
    // ! get the body geometry and wet/dry fraction
    hb.hh = hb.hhigh - hb.hlow;
    if(hb.hh <= 0.0) {
        // call errorhdr("hbeval: hhigh < hlow// !", "cont")
        printf("hbeval: hhigh < hlow// !, cont \n");
    }
    if(hb.tacr.n <= 0) {
        hb.acr = hb.vol/hb.hh; // ! assuming conservation of volume
    }
    else{
        hb.vol = hb.acr*hb.hh; // ! volume and mass changes // !
    }
    // end if
    hb.fwet = min(max(hb.clp.hlev-hb.hlow,0.0)/hb.hh, 1.0);
    if(hb.vol <= 0.0 || hb.acr <= 0.0) {
    // !     call errorhdr("hbeval: zero acr or vol// !","cont")
        hb.vol = 0.0; hb.acr = 0.0; hb.asf = 0.0; hb.ene = 0.0;
        hb.twet = hb.clp.tl; hb.twsf = hb.clp.tl;
        hb.tdry = hb.clp.tg; hb.tdsf = hb.clp.tg;
        hb.htmodw = "-"; hb.htmoddw = "-";
        hb.qgen = 0.0; hb.qtr = 0.0; hb.qrlx = 0.0;
        hb.qflxwet = 0.0; hb.qflxdry = 0.0;
        return;
    }

    // ! do the asf depending on ll
    if(hb.iasf == 0) {
        hb.asf = hb.vol*(2.0/hb.ll + 2.0/hb.hh);
    }
    

    // ! do the heat transfer
    if(hb.ihtr == 1) {
        // ! known coolant bulk state and properties:
        // !   ptot,ng,gname(),pg(),tg,tl,rhog(mixture),rhol,cpl,betal

        hwet = 0.0; hdry = 0.0; // ! initialize heat transfer coeff.
        kb = 2.0*hb.lam/hb.ftdrp/hb.ll; // ! gradient of qflux vs (tsf - tav) curve 

        // ! wet zone heat transfer
        hb.htmodw = '-'; is1 = 0;
        if(hb.fwet >0.0) {
            htcwet(time,dt,kb,hb,is1,cmode,hwet);
            hb.htmodw = cmode;
        }

        // ! dry zone heat transfer
        hb.htmoddw = '-'; is2 = 0;
        if(hb.fwet <1.0) {
            if(hb.icdn > 0 && hb.tdsf < hb.clp.tg) {
                htcdrycdn(time,dt,kb,hb,is2,hdry);
            }
            else {
                htcdry(time,dt,kb,hb,is2,hdry);
            }
            // ! if water injection is on and water mass>0, add wet htc assuming spray
            if(hb.fhtdrywet > 0.0 && (hb.clp.wvol > 0.0 && hb.clp.wmdi > 0.0)) {
                htcwet(time,dt,kb,hb,is2,cmode,hdrywet);
                hb.htmoddw = cmode;
                hdry = hdry*(1.0-hb.fhtdrywet) + hdrywet*hb.fhtdrywet;
            }
        }
    }
    else {// ! given heat transfer coeff. by input
        hwet = hb.htcwet;
        hdry = hb.htcdry;
    }

    // ! exchange between the dry and wet zones
    lrlx = min(hb.fwet,1.0-hb.fwet)*hb.hh*0.5;
    lrlx = max(lrlx, sqrt(hb.lam/hb.rho/hb.cp*100.0)*5.0);
    lrlx = min(lrlx, 0.5*hb.hh);
    // ! qrlx = frlx*acr*lam/lrlx*(tdry-twet) => included in the energy
    // ! equation with implicit scheme for temperatures. 

    // !  // ! give an arbitrary factor for heat transfer coefficients
    // !  hwet = hwet*hb.fhtr
    // !  hdry = hdry*hb.fhtr
    // !  modification is included in htcwet/dry

    // ! get tdry and twet by solving 2 linear equations
    // ! (solve the energy equations for tdry and twet)
    // ! a11*tdry + a12*twet = a13
    // ! a21*tdry + a22*twet = a23
    is3 = 0;
    if(hb.fwet > 0.0 && hb.fwet < 1.0) {
        x1 = hb.cp*hb.rho*hb.vol/dt;
        x2 = hb.frlx*hb.acr*hb.lam/lrlx;
        x3 = hb.ftdrp*hb.ll/2.0/hb.lam;
        a[0][0] = x1 + x2/(1.0-hb.fwet) + hdry*hb.asf/(1.0+x3*hdry);
        a[0][1] = -x2/(1.0-hb.fwet);
        a[1][0] = -x2/hb.fwet;
        a[1][1] = x1 + x2/hb.fwet + hwet*hb.asf/(1.0+x3*hwet);
        a[0][2] = hb.qgen + x1*hb.tdry + hdry*hb.asf/(1.0+x3*hdry)*hb.clp.tg ;
        a[1][2] = hb.qgen + x1*hb.twet + hwet*hb.asf/(1.0+x3*hwet)*hb.clp.tl ;
        gauss (a,2, is3 , neps, 1);
        hb.tdry = a[0][2]; hb.twet = a[1][2];
        // ! heat rate and fluxes for reference
        hb.tdsf = (hb.tdry + x3*hdry*hb.clp.tg)/(1.0+x3*hdry);
        hb.qflxdry = hdry*(hb.tdsf - hb.clp.tg);
        hb.twsf = (hb.twet + x3*hwet*hb.clp.tl)/(1.0+x3*hwet);
        hb.qflxwet = hwet*(hb.twsf - hb.clp.tl);
        hb.qrlx = x2*(hb.tdry-hb.twet);
    }
    else if (hb.fwet == 0.0) { // ! all dry
        x1 = hb.cp*hb.rho*hb.vol/dt;
        x3 = hb.ftdrp*hb.ll/2.0/hb.lam;
        hb.tdry = (hb.qgen + x1*hb.tdry + hdry*hb.asf/(1.0+x3*hdry)*hb.clp.tg) 
            /(x1 + hdry*hb.asf/(1.0+hdry*x3));
        hb.tdsf = (hb.tdry + x3*hdry*hb.clp.tg)/(1.0+x3*hdry);
        hb.qflxdry = hdry*(hb.tdsf - hb.clp.tg);
        hb.twet = hb.tdry; hb.twsf = hb.tdsf;
        hb.qflxwet = 0.0; hb.qrlx = 0.0;
    }
    else { // ! all wet
        x1 = hb.cp*hb.rho*hb.vol/dt;
        x3 = hb.ftdrp*hb.ll/2.0/hb.lam;
        hb.twet = (hb.qgen + x1*hb.twet + hwet*hb.asf/(1.0+x3*hwet)*hb.clp.tl) 
            /(x1 + hwet*hb.asf/(1.0+hwet*x3));
        hb.twsf = (hb.twet + x3*hwet*hb.clp.tl)/(1.0+x3*hwet);
        hb.qflxwet = hwet*(hb.twsf - hb.clp.tl);
        hb.tdry = hb.twet; hb.tdsf = hb.twsf;
        hb.qflxdry = 0.0; hb.qrlx = 0.0;
    }
    
    if((is1+is2+is3) > 0) {ist = 1;}

    // ! total heat transfer to the coolant (W), and the saturation temp. for reference
    hb.qtr = hb.qflxwet*hb.asf*hb.fwet + hb.qflxdry*hb.asf*(1.0-hb.fwet);
    hb.clp.tspt = tsatwaterf(hb.clp.ptot);
    // ! total energy of the body (base=tini) and the integrated heat generation (J)
    hb.ene = hb.rho*hb.acr*hb.hh*hb.cp 
        *((1.0-hb.fwet)*(hb.tdry-temp0)+hb.fwet*(hb.twet-temp0));
    hb.intqgen = hb.intqgen + hb.qgen*dt;
}
/*###############################################################################*/ 




































/*###############################################################################*/ 
// int main()
// {
//     htbody ht11;
//     ht11.htmodw = "C";
//     printf("TEST\n");
//     printf("TEST\n");
//     return 0;
// }
/*###############################################################################*/ 



#endif