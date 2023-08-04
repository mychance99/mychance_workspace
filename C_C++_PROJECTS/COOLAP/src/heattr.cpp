#ifndef HEATTR_H
#define HEATTR_H
#include <iostream>
// #include <vector>
#include <string>
#include <math.h>

using namespace std;

#include "consts.cpp"

/*###############################################################################*/ 
// ! A simple set of heat transfer correlations
// ! (natural convection, saturation pool boiling heat transfer, radiation)
// ! calculate heatflux
// !
// ! gas phase heat transfer
// ! htransfg (tsf, epsrad, tg, rhogf,cpgf,mugf,lamgf,betagf, 
// !      qflux, qcon, qrad)
// !
// ! under water heat transfer
// ! htransfl (tsf, epsrad, p, tsat, hfg,
// !     tl,rhol,cpl,mul,laml,sig, betal, 
// !     rhov,rhovf,cpvf,muvf,lamvf,betavf, 
// !     qflux, qboil, qcon, qrad, creg,
// !     tbi, tchf, tmfb, qchf, qmfb)
// !
// ! correlations for every heat transfer mode
// ! htnatconv (tsf,tc,rho,cp,lam,mu,beta,ll,grv, qcon)
// ! htrad (tsf,epsrad,tc, qrad)
// ! htnbkut (tsf, pres,tsat,hfg,rhov,rhol,cpl,laml,mul,sig, 
// !      qnb,tchf,qchf,iflg)
// ! htfbbere (tsf, tsat, hfg, rhol,sig, rhovf,cpvf,muvf,lamvf, 
// !      qfb, tminfb, qminfb, iflg)
/*###############################################################################*/ 

/*###############################################################################*/ 
void htrad (double tsf,double epsrad,double tc,double &qrad)
// ! radiation 
// ! stefan-boltzmann radiation law
// !
// ! addition of filmboiling and radiation HT
// ! sould be done by [BRO53]
// !    qtot = qfb + J*qrad, J=7/8
// !
// ! water emmisivity is assumed to be 1. because its near 1.
{
    qrad = epsrad*sigsb*max(pow(tsf,4)-pow(tc,4), 0.0);
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void htfbbere (double tsf, double tsat, double hfg, double rhol,double sig, 
               double rhovf,double cpvf,double muvf,double lamvf, 
               double &qfb, double &tminfb, double &qminfb, int &iflg)
// ! film-boiling heat transfer model 
// ! by berenson [BER61]
// ! "min. film boiling temperature is also given.
// !  and radiative heat transf. for horizonal upward surface
// !
// ! behavior of this subroutine
// !   calc. qminfb and tminfb
// !   if tsf < tminfb
// !       qfb=qminfb => give back (tminfb, qfb(tminfb)) for
// !       iflg = 1    connection with transition boil
// !   else
// !       calc. qfb at tsf =>  in film boiling regime
// !       iflg = 0
// !   end if
// !
// !   radiation effect should be externally considered
{
// ! definition of nusselt num. and heat transf. coeff.
// !   nufb = hfb*ll/lamvf  (ll: laplace length)
// !   hfb = qfb/(tsf-tsat)  
// !    => qfb = hfb*(tsf-tsat) = nufb*lamvf/ll*(tsf-tsat)
// !  (tsat may be replaced by tl that is not pressure dependent) 
// ! definition of min. film boiling temp.
// !   tminfb = minimum tsf at which film boiling can sustain
// !
// !   normal liq. properties are used here for brevity  
// !   although precisely they are film properties at (tsat+tl)/2
    double prvf, hfg1, dtminfb, 
       sp1, ll, ccmin, ccht, gr1, nu, 
       fnew, fdnew, delnew, xx, yy, zz, tw;
    int inew, i;
    prvf = cpvf*muvf/lamvf;
    // ! ========================================================
    // ! Berenson min. film boiling temp. correlation [BER61]
    // ! ========================================================
    // !  original correlation
    // !  qminfb = C*hfg*rhovf
    // !       /(rhovf**2/sig/grav/(rhol-rhovf))**0.25
    // !       /((rhol-rhovf)/rhovf)**0.5
    // !  C=0.09   : const.
    // !  get tminfb along qfb vs tsf curve by Berenson
    // !
    ccmin = 0.09;
    qminfb = ccmin*hfg*rhovf 
        /pow(max( pow(rhovf,2)/max(sig,1.e-4)/grav 
                /max(rhol-rhovf,10.0), 1e-12 ),0.25)
        /pow(max( (rhol-rhovf)/max(rhovf,1e-6) , 1e-12 ),0.5);
    // !
    ccht = 0.425;
    ll = sqrt( max( sig/grav/max(rhol-rhovf, 10.0), 0.0) ); // ! laplace length
    gr1 = grav*rhovf*max(rhol-rhovf, 10.0)*pow(ll,3)  // ! mod. grashof no.
        /max( pow(muvf,2), 1.e-16 );
    xx = pow(( qminfb*ll/ccht/max(lamvf, 1.e-12) ),4);
    yy = gr1*prvf*hfg/max(cpvf, 1e-6);
    zz = 0.5*gr1*prvf;
    dtminfb = 500.0;
    inew = 0;

    for (i=1; i<=20; i++) // do i=1,20
    {
        fnew = zz*pow(dtminfb,4) + yy*pow(dtminfb,3) - xx;
        fdnew = 4.0*zz*pow(dtminfb,3) + 3.0*yy*pow(dtminfb,2);
        if(fabs(fdnew) < 1.e-12) {
            // write(*,'(a,1pe12.4)')  &
            //     'htfbbere: newton meth. singlar // ! dtminfb=', dtminfb
            printf("htfbbere: newton meth. singlar // ! dtminfb= %f", dtminfb);
            delnew = 500.0;
            if(inew != 0) {exit(EXIT_FAILURE);}// stop
            inew = 1;
        }
        else {
            delnew = -fnew/fdnew;
            dtminfb = dtminfb + delnew;
        }
        if(fabs(delnew) < 1.e-6) {exit(EXIT_FAILURE);}// exit
    }
    if(fabs(delnew) > 1.0) {
        // write(*,'(a,1p,2(e12.4,1x))')  &
        //     'htfbbere: dtminfb not converged dtminfb,delnew=', &
        //     dtminfb, delnew
        printf("htfbbere: dtminfb not converged dtminfb,delnew=%f %f", dtminfb, delnew);
    }
    tminfb = tsat + max(dtminfb, 0.0);
    // !
    // ! ========================================
    // ! switch heat source temperature and flag
    // ! ========================================
    if( tsf < tminfb ) {
        iflg = 1;     // ! to calculate film boiling t,q" for interpolation to 
        tw = tminfb;  // ! the transition boiling q"
    }
    else {
        iflg = 0;
        tw = tsf;
    }
    // !
    // ! =============================================================
    // ! Berenson film boiling correlation for 
    // ! horizontal upward surface [BER61]
    // ! =============================================================
    // !
    hfg1 = hfg + 0.5*cpvf*max(tw-tsat,0.0);
    sp1 = max( cpvf*max(tw-tsat,1.0)/hfg1/prvf, 1.e-6 );
    nu = ccht*pow(max( gr1/sp1, 0.0 ),0.25);

    // ! heat flux by conductive/convective component of film boiling
    qfb = nu*lamvf/max(ll, 1.e-4)*max(tw-tsat,0.0);
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void htnbkut (double tsf, double pres,double tsat,double hfg,double rhov,
              double rhol,double cpl,double laml,double mul,double sig,
              double &qnb,double &tchf,double &qchf,int &iflg)
// ! nucleate boiling model
// ! kutateladze [KUT52]
// !
// ! behavior of this subroutine
// !   calc. qchf [ZUB61] and Tchf
// !   if tsf > tchf
// !       set qnb=qchf => give back (tchf,qchf) for connection with 
// !       iflg = 1        transition boil
// !   else
// !       calc. qnb
// !       iflg = 0
// !   end if
{
    double prl, aa, bb, cc;
    // ! CHF heat flux [ZUB61]
    qchf = 0.131*hfg*rhov/pow(max( pow(rhov,2)/max(sig,1e-4)/grav/ 
        max(rhol-rhov, 10.0), 1e-12 ),0.25);
    // ! CHF temperature based on qchf and Kutateladze correlation
    aa = max( sig/grav/max(rhol-rhov,10.0), 0.0 );
    bb = max(pres*rhol/max(sig*hfg*rhov*mul, 1e-6), 0.0);
    cc = max(7e-4 * laml * pow(prl,0.35), 0.0);
    tchf = tsat + pow(qchf,0.3)/max( cc * pow(aa,0.2) * pow(bb,0.7), 1e-6 );

    if( tsf > tchf ) {
        iflg = 1;
        qnb = qchf;
    }
    else {
        iflg = 0;
        // ! Kutateladze boiling curve
        qnb = pow(cc,3.33) * pow(max(tsf-tsat, 0.0),3.33)
           * pow(aa,0.667) * pow(bb,2.33);
        qnb = min(qnb, qchf);
    }
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void htnatconv (double tsf,double tc,double rho,double cp,double lam,
                double mu,double beta,
                double ll,double grv,
                double &qcon)
// ! #################################################################
// ! natural convection (Churchil-Chu) (Shoji text book)
// !  if gravity is omitted, 9.8 is used
// !  if length scale is omitted, 1 is used (the result does not
// !  significantly depends.)
// ! #################################################################
{
    // ! tsf : surface temp.(K)      tc: coolant temp. (K)
    // ! ll: heating surface length (basically vertical) (m)
    // ! rho: coolant density(kg/m3)  cp: specific heat(J/kgK)
    // ! lam: heat conductiviety(W/mK)  mu: viscosity(Pa)
    // ! beta: volume expansion coefficient (1/K) =(1/v0)(dv/dt)=(-1/rho0)(drho/dt)
    // ! recommend the usage of film properties at (tsf+tc)/2 especially for gas
    // ! with which the temperature difference between the surface and the bulk 
    // ! fluid can be large
    double g, l, pr, ra, nu;

    // if(grv != 0) {
        g = grv;
    // }
    // else {
        // g = grav;
    // }
    if(ll != 0) {
        l = ll;
    }
    else {
        l = 1.0; // ! dependence on ll is not so significant, so...
    }
    pr = (mu/rho)/(lam/cp/rho);
    ra = beta*g*fabs(tsf-tc)*pow(ll,3)/pow((mu/rho),2)*pr;
    nu = pow((0.825 + 0.387*pow(ra,(1.0/6.0))/pow((1.0+pow((0.492/pr),(9.0/16.0))),(9.0/27.0))),2);
    qcon = nu*lam/ll*(tsf - tc);
}           
/*###############################################################################*/ 


/*###############################################################################*/ 
void htransfl (double tsf, double epsrad, double p, double tsat, double hfg,
               double tl,double rhol,double cpl, double mul, double laml, 
               double sig, double  betal, double  rhov, double rhovf, double cpvf, 
               double muvf, double lamvf, double betavf, 
               double &qflux, double &qboil, double &qcon, double &qrad, string &creg, 
               double &tbi, double &tchf, double &tmfb, double &qchf, double &qmfb)
{
    double qfb, qnb, qchf0=0.0,qmfb0=0.0,tbi0=0.0,tmfb0=0.0,tchf0=0.0, ll = 1.0;
    int ifb, inb;

    // ! radiation HT
    htrad (tsf,epsrad, tl, qrad);

    // ! boiling inception temp.
    // ! simply tsf>tsat is used because surface temp. is given.
    tbi0 = tsat;

    ifb = -1;  // ! initialize flags
    inb = -1;

    if (tsf < tbi0) {  // ! low temp. no-boiling regime
        qboil = 0.0;
        creg = "c";
    }
    else{                  // ! now, boiling is considered
        htfbbere           // ! film boiling model
            (tsf, tsat,hfg, rhol,sig, rhovf,cpvf,muvf,lamvf, 
            qfb, tmfb0, qmfb0,ifb);
        htnbkut            // ! nuc. boiling model
            (tsf, p,tsat,hfg,rhov,rhol,cpl,laml,mul,sig, 
            qnb, tchf0, qchf0, inb);
        if (ifb == 0) {    // ! film boiling
            qboil = qfb;
            creg = "f";
        }
        else {                // ! not film boiling
            if(inb == 0) {    // ! nuc. boil.
                qboil = qnb;
                creg = "n";
            }
            else {            // ! not nuc. boil. it"s transition
                // ! transition boiling : linear interplation between
                // ! (tchf,qchf)--(tminfb,qmfb)
                if(tmfb0 > tchf0) {
                    qboil = qchf0-(qchf0-qmfb0)/(tmfb0-tchf0)*(tsf-tchf0);
                }
                else {       // ! this should not happen, but...
                    qboil = 0.5*(qmfb0+qchf0);
                }
            creg = "t";
            }
        }
    }

    // ! convection
    htnatconv(tsf,tl,rhol,cpl,laml,mul,betal,ll,grav, qcon);

    // ! synthesis of the heat transfer (continuity of qflux is required)
    if(creg == "c") {  // ! convection 
        qflux = qcon + qrad;
    }
    else if(creg == "n") {    // ! nucleate boiling
    // !     qflux = max(qboil, qcon) + 7.0/8.0*qrad 
        qflux = max(qboil,qcon) + qrad;
    }
    else {
    // !     qflux = qboil + 7.0/8.0*qrad 
        qflux = qboil + qrad;
    }


    // if present
    tbi = tbi0;
    qmfb = qmfb0;
    qchf = qchf0;
    tchf = tchf0;
    tmfb = tmfb0;

}
/*###############################################################################*/ 

#endif