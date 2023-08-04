#ifndef DBRCINTEG_H
#define DBRCINTEG_H

using namespace std;
#include <math.h>

#include "dbrcool.cpp"
#include "hotcan.cpp"
#include "wrsttab.cpp"

/*###############################################################################*/ 

/*###############################################################################*/ 
void ifdbrccan(meltdebr &md,htcan &cn,double time, double dt)
{
    double qrate;

    readmhex(md.hx,time,dt, qrate); // ! qrate for time..time+dt average (W)

    cn.qsrcex = cn.qsrcex + qrate;

}
/*###############################################################################*/ 

/*###############################################################################*/ 
void initdbrc(meltdebr &md)
{
    // ! initialization of meltdebr object
    // ! this should be called after initialization of hotcan where meltdebr is put
    // ! and reading the can's water pool geometry info

    // ! geometric condition input check
    // ! water pool area and radius

    if(md.plg.csar > 0.0) {
        md.plg.rr = sqrt(md.plg.csar/pi);
    }
    else if(md.plg.rr > 0.0) {
        md.plg.csar = pi*pow(md.plg.rr,2);
    }
    else
    {
        // call errorhdr(mess="initdbrc: need md.plg.csar or md.plg.rr",act="stop")
        printf("initdbrc: need md.plg.csar or md.plg.rr\n");
        exit(EXIT_FAILURE);
    }
    // ! water height
    if(md.plg.hhp > 0.0) {
        if(md.plg.hhp > md.plg.hh0) {
        // call errorhdr(mess="initdbrc: md.plg.hhp>md.plg.hh0; give hhp0 if pool level>hh0",act="stop")
        printf("initdbrc: md.plg.hhp>md.plg.hh0; give hhp0 if pool level>hh0\n");
        exit(EXIT_FAILURE);
        }
        else if (md.plg.hhp0 > 0.0) {
        // !  call errorhdr(mess="initdbrc: md.plg.hhp0 overrides md.plg.hhp",act="cont")
        md.plg.hhp = min(md.plg.hhp0,md.plg.hh0);
        }
        else
        {
        md.plg.hhp0 = md.plg.hhp;
        }
    }
    else if(md.plg.hhp0 > 0.0) {  // ! this is the normal input style
        md.plg.hhp = min(md.plg.hhp0,md.plg.hh0);
    }
    else {
        // call errorhdr(mess="initdbrc: need pool height md.plg.hhp0",act="stop")
        printf("initdbrc: need pool height md.plg.hhp0\n");
        exit(EXIT_FAILURE);
    }
    // ! melt release position
    if(md.plg.hh0 < 0.0){
        // call errorhdr(mess="initdbrc: need melt release position md.plg.hh0",act="stop")
        printf("initdbrc: need melt release position md.plg.hh0\n");
        exit(EXIT_FAILURE);
    }
    // ! check total pressure (only variable needed for gas phase)
    if(md.clp.ptot <= 0.0) {
        // call errorhdr(mess="initdbrc: pressure not given",act="stop")
        printf("initdbrc: pressure not given\n");
        exit(EXIT_FAILURE);
    }

    initmprop(md.mppb);
    initdbed(md.db,md.plg.rr);
    initmlump(md.ml);
    initmhex(md.hx);

    // ! initialize clp (coolant prop data) => not needed actually, but for future
    md.clp.ng = 3; // ! water vapor, non-condensible gas (air) and hydrogen
                    // ! due to the spec. of hotcan (3 species)
                    // ! => EXTENTION of hotcan NEEDED FOR MORE SPECIES (e.g. CO, CO2..)
    md.clp.gname[1] = "h2o";
    md.clp.gname[2] = "air";
    md.clp.gname[3] = "h2";



}
/*###############################################################################*/ 


/*###############################################################################*/ 
void ifcandbrc(meltdebr &md, htcan cn,int init=0)
{
    int i1, ini;
    double x0,x1,x2,x3;

    ini = 0;
    // if(present(init)) ini = init
    ini = init; 

    // ! coolant properties: only ptot and relevant water/vapor properties are needed
    md.clp.ptot = cn.ptot;
    md.clp.tl=cn.tl; md.clp.tg=cn.ts;
    md.clp.rhol=cn.rhol; 
    md.clp.cpl=cn.cpl; 
    md.clp.betal=cn.betal;
    md.clp.hlev=cn.hl;
    md.clp.wvol=cn.voll; md.clp.wmdi=cn.mdil;

    // ! set pool geometric vars
    md.plg.hhp0 = md.clp.hlev;  // ! water pool depth
    if(md.plg.hhp0 > 0.0 && (md.plg.csar <= 0.0 && md.plg.rr <= 0.0)) {
        md.plg.csar = md.clp.wvol/md.plg.hhp0;
    }
    md.plg.hhp = min(md.plg.hhp0,md.plg.hh0);
    md.plg.hhf = md.plg.hh0 - md.plg.hhp;

    // ! calculate coolant physical properties not obtained from hotcan
    md.clp.tspt = tsatwaterf(md.clp.ptot);
    md.clp.hfg = hfgwaterf(md.clp.ptot);
    md.clp.tg = md.clp.tspt;  // ! gas=steam: only consider inside water pool
    // i1 = wrsteamtab(md.clp.ptot,md.clp.tspt,1,md.clp.rhov,x0,x1,x2,x3,md.clp.cpv,i1);
    wrsteamtab(md.clp.ptot,md.clp.tspt,1,md.clp.rhov,x0,x1,x2,x3,md.clp.cpv,i1);
    md.clp.lamv = tconwaterf(md.clp.rhov,md.clp.tspt);
    md.clp.muv = viscwaterf(md.clp.rhov,md.clp.tspt);
    md.clp.sig = stenwaterf(md.clp.tspt);
    md.clp.laml = tconwaterf(md.clp.rhol,md.clp.tl);
    md.clp.mul = viscwaterf(md.clp.rhol,md.clp.tl);

    // ! calculate water pool macroscopic balance
    md.wpoolmass = md.clp.rhol * md.clp.wvol;
    md.wpoolenth = md.wpoolmass * (cn.el + md.clp.ptot/md.clp.rhol);
    x0 = (cn.ev + cn.pv/cn.rhov)* cn.mv;
    if(ini == 1) {
        md.wpoolmass0 = md.wpoolmass;
        md.enthvap0 = x0;
    }
    md.wpoolmassevap = md.wpoolmass0 - md.wpoolmass;
    md.wpooleneevap = x0 - md.enthvap0;
}
/*###############################################################################*/ 




// !Calculating jet initial velcoity
void initialjetvel(meltdebr &md, double time, double dt)
{

  double volhemisphere, volcorium, hvessel, re, ff, pdif, timei, vv;

  // !  md.plg.hhf = md.plg.hh0 - md.plg.hhp  // !Corium free fall height
  // !  ttabeval(md.djin, time, md.mj.d0)  // !Corium jet diameter from time table
  // !  ttabeval(md.tjin, time, md.mj.temp0) // !Corium jet temperature from time table

  // !  md.mj.tempi = md.mj.temp0 // !Corium initial temperature
  // !  md.mj.mpp.tm = md.mj.tempi // !Corium temperature for calculating other properties
  // !  getmprop(md.mppb, md.mj.mpp)
  if(md.mj.pin>md.clp.ptot) {
    pdif=md.mj.pin-md.clp.ptot; // !Pressure difference between in-vessel and containment
  }
  else {// !It must be P_in > P_out
    md.mj.pin=md.clp.ptot;
    pdif=0;
  }
  // endif
  volhemisphere = 2.0/3.0*pi*pow(md.mj.r0,3)*abs(cos(abs(asin(md.mj.d0/2.0/md.mj.r0)))); // !Volume of vessel lower head
  volcorium = md.mj.mretot/md.mj.mpp.rhom; // !Volume of corium in vessel 


  // !jmod 1=using given velocity, 2=bernoulli without considering friction(global failure), 3=bernoulli with considering friction(penetration failure)
  if(md.mj.jmod==1) {
    ttabeval(md.vjin, 0.0, vv);
    ttabeval(md.vjin, time, md.mj.v0);
    timei = -vv/grav + sqrt(pow(vv,2)+2.0*grav*(md.plg.hh0-md.plg.hhp0))/grav; // ! time of water surface contact
  }
  if(time>timei) {
    ttabeval(md.vjin, time-timei, md.mj.v0);
  }
  // endif
  else if(md.mj.jmod==2 || md.mj.v0==0) {
    if(volhemisphere<volcorium) {
      hvessel = sqrt(pow(md.mj.r0,2)-pow((md.mj.d0/2.0),2)) +(volcorium-volhemisphere)/pi/pow(md.mj.r0,2);
      md.mj.v0 =sqrt(2.0*pdif/md.mj.mpp.rhom+2.0*grav*hvessel);
      md.mj.mretot = md.mj.mretot - pi/4*pow(md.mj.d0,2)*md.mj.mpp.rhom*md.mj.v0*dt;
    }
    else if(volcorium<=0.0) {
      hvessel = 0.0;
      md.mj.v0 =0.0;
      md.mj.mretot = 0.0;
    }
    else {
      hvessel = sqrt(pow(md.mj.r0,2)-pow((md.mj.d0/2.0),2))-md.mj.r0*(cos(asin(md.mj.d0/2.0/md.mj.r0))-3.0*volcorium/2.0/pi/pow(md.mj.r0,3));
      md.mj.v0 =sqrt(2.0*pdif/md.mj.mpp.rhom+2.0*grav*hvessel);
      if(md.mj.llj==0) {
        md.mj.mretot = md.mj.mretot; // !mass err correct
      }
      else {
        md.mj.mretot = md.mj.mretot - pi/4*pow(md.mj.d0,2)*md.mj.mpp.rhom*md.mj.v0*dt;
      }
    }
  }
  else if(md.mj.jmod==3) {
    re=(md.mj.mpp.rhom*md.mj.v0*md.mj.d0)/md.mj.mpp.mum;
    if(re<10e+4) {
      ff=0.079/pow(re,0.25);                       // !fanning friction factor (Balsius correlation)
    }
    else if(re>=10e+4 && re<10e+5) {
      ff=(0.0014+0.125/pow(re,0.32)-0.079/pow(re,0.25))/(pow(10,5)-pow(10,4))*(re-pow(10,4))+0.079/pow(re,0.25); // !interpolation of two fff models
    }
    else {
      ff=0.0014+0.125/pow(re,0.32);                  // !fanning friction factor (koo correlation)
    }
    if(volhemisphere<volcorium) {
      hvessel = sqrt(pow(md.mj.r0,2)-pow((md.mj.d0/2.0),2)) +(volcorium-volhemisphere)/pi/pow(md.mj.r0,2);
      md.mj.v0 =sqrt(2.0*pdif/md.mj.mpp.rhom+2.0*grav*(hvessel+md.mj.lb))/sqrt(4*ff*md.mj.lb/md.mj.d0+1+md.mj.kk);
      if(md.mj.llj==0) {
        md.mj.mretot = md.mj.mretot; // !mass err correct
      }
      else {
        md.mj.mretot = md.mj.mretot - pi/4*pow(md.mj.d0,2)*md.mj.mpp.rhom*md.mj.v0*dt;
      }
    }
    else if(md.mj.mretot<=0.0) {
      hvessel = 0.0;
      md.mj.v0 =0.0;
      md.mj.mretot = 0.0;
    }
    else {
      hvessel = sqrt(pow(md.mj.r0,2)-pow((md.mj.d0/2.0),2))-md.mj.r0*(cos(asin(md.mj.d0/2.0/md.mj.r0))-3.0*volcorium/2.0/pi/pow(md.mj.r0,3));
      md.mj.v0 =sqrt(2.0*pdif/md.mj.mpp.rhom+2.0*grav*(hvessel+md.mj.lb))/sqrt(4*ff*md.mj.lb/md.mj.d0+1+md.mj.kk);
      if(md.mj.mretot<=pi/4*pow(md.mj.d0,2)*md.mj.mpp.rhom*md.mj.v0*dt) {
        md.mj.v0 = md.mj.mretot/(pi/4*pow(md.mj.d0,2)*md.mj.mpp.rhom*dt);
      }
      if(md.mj.llj==0) {
        md.mj.mretot = md.mj.mretot; // !mass err correct
      }
      else
        md.mj.mretot = md.mj.mretot - pi/4*pow(md.mj.d0,2)*md.mj.mpp.rhom*md.mj.v0*dt;
    }
  }
}
















/*###############################################################################*/ 
void dbrcevol(meltdebr &md, double time, double dt) // << by juwook 2022
// ! time evolusion for one time step
{
    md.plg.hhf = md.plg.hh0 - md.plg.hhp;
    ttabeval(md.djin, time, md.mj.d0);
    ttabeval(md.vjin, time, md.mj.v0);
    if (time == 0) {
        ttabeval(md.tjin, time, md.mj.temp0);
    }
    md.mj.tempi = md.mj.temp0; // corium initial temperature
    md.mj.mpp.tm = md.mj.tempi; // corium temp. for calculating other properties
    getmprop(md.mppb, md.mj.mpp);

    // calc. init. vel. v0 by pressure difference
    // initialjetvel(md, time, dt);


    // ! melt jet evolution
    // ! mass/energy deposition into other components are called in it
    mjetevol(md,time,dt);

    // ! evolution in debris bed and lump
    // ! heat release and energy conservation; remelt and merge to lump
    dbrlmpevol(md,time,dt);

}
/*###############################################################################*/ 





#endif