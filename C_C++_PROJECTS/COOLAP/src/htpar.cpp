#ifndef HTPAR_H
#define HTPAR_H

/*###############################################################################*/ 
#include <iostream>
#include <math.h>
using namespace std;

// #include "consts.cpp"
/*###############################################################################*/ 

/*###############################################################################*/ 
void htfboilpar (double pres, double alp, double dd, double vl, double tl, 
                 double tsf, double epsrad, double tsat, double hfg, double rhol, 
                 double cpl, double mul, double laml, double sig, double rhovf, 
                 double cpvf, double lamvf, double muvf, 
                 double &qfb, double &qrad, double &tminfb, int &iflg)
{
    /* ! film-boiling heat transfer model
    ! by liu & theofanous [LIU95a]
    ! "real" version with Kondo [KON95] min. film boiling temperature
    !  and radiative heat transf. for particle model
    !
    ! behavior of this subroutine
    !   calc. qrad
    !   calc. tminfb
    !   if tsf < tminfb
    !       calc. qfb at tminfb  => give back (tminfb, qfb(tminfb)) for
    !       iflg = 1                 connection with transition boil
    !   else
    !       calc. qfb at tsf =>  in film boiling regime
    !       iflg = 0
    !   end if
    !
    ! definition of nusselt num. and heat transf. coeff.
    !   nufb = hfb*dp/lamvf
    !   hfb = qfb/(tsf-tl)  
    !   # precisely it is (tsf-tsat) but tl was taken for stability
    !    => qfb = hfb*(tsf-tl) = nufb*lamvf/dp*(tsf-tl)
    !
    ! definition of min. film boiling temp.
    !   tminfb = minimum tsf at which film boiling can sustain
    !
    !  notes:
    !   addition of filmboiling and radiation HT (lator)
    !   sould be done by [BRO53]
    !      qtot = qfb + J*qrad, J=7/8
    !      qrad = eps*sigb*(tsf^4 - tl^4)
    !
    !   normal liq. properties are used here for brevity  
    !   although precisely they are film properties at (tsat+tl)/2 */

    /* real(kind(1d0)),intent(in) :: &
        pres,    & ! ambient pressure 
        alp,     & ! void frac.
        dd,      & ! dia. of hot sphere
        vl,      & ! liquid-hotsphere relative velocity
        tl, tsf, & ! temperatures (liq., hot surface)
        epsrad,  & ! melt emissivity
        tsat,hfg,& ! sat. temp. and latent heat at pres
        sig, cpl, rhol, mul, laml, & ! liq. props.
        rhovf, cpvf, muvf, lamvf     ! vapor film props. at tvf=(tsat+tsf)/2
    real(kind(1d0)),intent(out) :: &
        qfb,    & ! heat flux by film boiling (Liu-Theo, no radiation)
        qrad,   & ! heat flux by radiation (Stefan-Boltzmann)
        tminfb    ! min film boiling temp. (Kondo)
    integer,intent(out) :: iflg
        ! =0 film boiling
        ! =1 not film boiling, give back qfb at tminfb
        !    for connection to transition boiling
    real(kind(1d0)) :: &
        hfg1, prl, prvf,  &
        sc1, sp1, d1, aar, kk, rr, cc, bb, aa, ee, mmc, &
        rel, fr, fff, kkc, ffa, ccnp, nup, nus, nuf, nu, nueps,  &
        fnew, fdnew, delnew, xx, yy, zz, &
        tw, fard, tmp1, tmp2, tmp3, tmp4
    real(kind(1d0)),parameter :: tcr=647.3d0
    integer :: itmp1, inew, i */
    double hfg1, prl, prvf, sc1, sp1, d1, aar, kk, rr, cc, bb, aa, ee, mmc, 
           rel, fr, fff, kkc, ffa, ccnp, nup, nus, nuf, nu, nueps, fnew, fdnew, 
           delnew, xx, yy, zz, tw, fard, tmp1, tmp2, tmp3, tmp4;
    double const tcr = 647.3;
    int itmp1, inew, i;

    prvf = cpvf*muvf/lamvf;
    prl = cpl*mul/laml;

    // ! ========================================================
    // ! stefan-boltzmann radiation law
    // ! ========================================================
    // ! water emmisivity is assumed to be 1. because its near 1.
    qrad = epsrad*sigsb*max(pow(tsf,4)-pow(tl,4), 0.);

    // ! consider void 
    // !   absourption into the local gas is dropped.
    // !   actually radiation go through gas and reach beyond cells
    // !   but HT beyond cells does not fits in the present framework.
    // !   if it is deposited on the local gas, it will jump the local
    // !   gas temp.  so, it is simply dropped. i.e. void reduces qrad.
    // !
    // ! alp=0-0.3:1-alp, alp=0.3-0.95: 0.7-0
    // !      fard = min(1.0-alp, max(0.7-1.07*(alp-0.3),0.))
    // ! alp<0.3:1, alp=0.3-0.95: 1->0
    // !      fard = min(1., max(0., 1.-(alp-0.3)/(0.95-0.3)))
    // !
    // ! alp<0.6:1, alp=0.6-0.95: 1->0
    fard = min(1., max(0., 1.-(alp-0.6)/(0.95-0.6)));
    qrad = qrad*fard;

    // ! ========================================================
    // ! Kondo et al. min. film boiling temp. correlation [KON95]
    // ! ========================================================
    // !  original correlation
    // !   tminfb = tsat + C*(27/32*tcr - tsat) +
    // !           (laml/lamv)*(tsat-tl)*Nuc/(dd/delmin+Nur)
    // !   C=0.6   : const.
    // !   delmin=0.1e-3 : min thickness of vapor film
    // !   Nuc : convection nusselt no. around vapor film surface (liq.side)
    // !   Nur : radiative nusselt no. around sphere
    // !  the equation is solved by newton method because Nur involves
    // !  tminfb**3
    xx = epsrad*sigsb*dd/lamvf;
    yy = dd/0.1e-3;
    zz = -tsat-0.6*max(27./32.*tcr-tsat, 0.);
    rel = rhol*vl*dd/mul;
    bb = laml/lamvf*max(tsat-tl,0.)* 
        min(2.+0.6*pow(rel,0.5)*pow(prl,0.333333), 100.);
    tminfb = 1000.;
    inew = 0;
    i=1;
    while (i<21){
        fnew = (xx*pow(tminfb,3)+yy)*(tminfb+zz)-bb;
        fdnew = 4.*xx*pow(tminfb,3)+3.*xx*zz*pow(tminfb,2)+yy;
        if(fabs(fdnew) < 1.e-12){
            printf("*htfboilpar* newton meth. singlar// ! tminfb=%f", tminfb);
            delnew = 500.;
            if(inew != 0) {break;}
            inew = 1;
        }
        else{
            delnew = -fnew/fdnew;
            tminfb = tminfb + delnew;
        }
        // end if
        if(fabs(delnew) < 1.e-6) {break;}
        i+=1;
    }
    // end do

    if(fabs(delnew) > 1.){
        printf ("'*htfboilpar* tminfb not converged// ! tminfb,delnew=%f, %f'",tminfb, delnew);
    }
    // end if
    tminfb = max(tminfb, tsat);

    // ! ========================================
    // ! switch heat source temperature and flag
    // ! ========================================
    if( tsf < tminfb ){
        iflg = 1;
        tw = tminfb;
    }
    else{
        iflg = 0;
        tw = tsf;
    }
    // end if

    // ! =============================================================
    // ! liu-theofanous film boiling correlation for a sphere [LIU95a]
    // ! =============================================================
    hfg1 = hfg + 0.5*cpvf*max(tw-tsat,0.);
    sc1 = cpl*max(tsat-tl,0.)/hfg1/prl;
    sp1 = cpvf*max(tw-tsat,1.)/hfg1/prvf;
    d1 = dd/sqrt( max( sig/(grav*(rhol-rhovf)), 1.e-12) );
    d1 = max(d1, 1.e-6);
    aar = grav*(rhol-rhovf)*pow(dd,3)/(pow(muvf,2)/rhovf);
    kk = rhol/max(rhovf, 1.e-12);
    rr = sqrt(muvf*rhovf/max(mul*rhol, 1.e-12));
    cc = 0.5*pow(rr,2)*sp1*prl;
    bb = -4./27.*pow(sc1,2) + 2./3.*sp1*prl*sc1 
        -32./27.*sp1*prl*pow(rr,2) + 0.25*pow(sp1,2)*pow(prl,2)
        +2./27.*pow(sc1,3)/pow(rr,2);
    bb = max(bb,0.);
    aa =  1./27.*pow(sc1,3) + 1./3.*pow(rr,2)*sp1*prl*sc1  
        +0.25*pow((rr*sp1*prl),2);
    ee = pow(max(aa+cc*pow(bb,0.5), 0.),0.33333)
        + pow(max(aa-cc*pow(bb,0.5), 0.),0.33333) 
        +1./3.*sc1;
    mmc = pow(ee,3)/(1.+ee/sp1/prl)/pow((rr*prl*sp1),2);

    rel = vl*dd*rhol/mul;
    rel = min(max(rel, 0.), 1.5e5);  // ! give limit from LIU95a data
    fr = pow(vl,2)/grav/max(dd, 1.e-6);
    fr = max(fr, 0.);
    fff = 1.-0.2/(1.+max(pow(fr,0.5)-1., 0.));

    if(d1 <= 0.14){
        kkc = 0.5/pow(d1,0.25);}
    else if(d1 <= 1.25){
        kkc = 0.86/(1.+0.28*d1);}
    else if(d1 <= 6.6){
        kkc = 2.4*d1/(1.+3.*d1);}
    else{
        kkc = 0.47*pow(d1,0.25);}
    // end if

    // ! pool : solve nup/(1+nup/2)=C -> nup = (c+sqrt(c^2+8c))/2
    ccnp = kkc*pow(max( aar/sp1*mmc, 0.),0.25);
    nup = 0.5*(ccnp + sqrt(max(pow(ccnp,2)+8.*ccnp,0.)));
    // ! conv. sat.
    nus = 0.5*pow(rel,0.5)*mul/muvf*pow(max(kk*pow(rr,4)/sp1,0.),0.25);
    // ! conv. sub.
    nuf = nus + 0.072*pow(rel,0.77)*pow(prl,0.5)*mul/muvf*sc1/sp1;
    // ! general nu
    nu = pow(max(pow(nup,5) + pow((fff*nuf),5), 0.),0.2);

    // ! two phase correction
    if(alp < 0.65){
        ffa = 1.;}
    else if(alp < 0.7){
        ffa = 1.-5.2*max(alp-0.65, 0.);}
    else if(alp < 0.93){
        ffa = pow((1.-alp),0.25);}
    else{
        ffa = 25.7*max(0.95-alp, 0.);}
    // end if

    nu = nu*ffa;

    // ! heat flux by conductive/convective component of film boiling
    qfb = nu*lamvf/max(dd, 1e-3)*max(tw-tl,0.);  // ! give limiter to dd


    return;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void htnucboil (double pres, double alp, double tl, double tsf, double tsat, 
                double hfg, double  rhov, double rhol, double cpl, double laml, 
                double mul, double sig, 
                double &qnb, double &tchf, int &iflg) 
{

    /* ! nucleate boiling model
    ! kutateladze [KUT52]
    !
    ! note: this subroutine canbe shared with melt pool model
    !
    ! behavior of this subroutine
    !   calc. qchf [ZUB61] and Tchf
    !   if tsf > tchf
    !       set qnb=qchf => give back (tchf,qchf) for connection with 
    !       iflg = 1        transition boil
    !   else
    !       calc. qnb
    !       iflg = 0
    !   end if */

    /* real(kind(1d0)),intent(in) :: &
        pres,    & ! pressure
        alp,     & ! void frac.
        tl, tsf, & ! temperatures
        tsat,hfg,& ! sat. temp. and latent heat at pres
        rhov, rhol, cpl, laml, mul, sig ! properties
    real(kind(1d0)),intent(out) :: &
        qnb,   & ! nuc. boiling heat flux (=qchf if chf condition)
        tchf     ! chf temperature
    integer,intent(out) :: iflg
        ! =0  nucleate boiling
        ! =1  exceed chf, give back qchf for connection to transition boiling
    real(kind(1d0)) :: prl, qchf, aa, bb, cc, fa */
    double prl, qchf, aa, bb, cc, fa;

    prl = cpl*mul/laml;

    // ! CHF heat flux [ZUB61]
    qchf = 0.131*hfg*rhov/pow(max( pow(rhov,2)/max(sig,1.e-4)/grav/ 
            max(rhol-rhov, 10.) , 1.e-12 ),0.25);
    // ! CHF temperature based on qchf and Kutateladze correlation
    aa = max(sig/grav/max(rhol-rhov,10.), 0.0);
    bb = max(pres*rhol/max(sig*hfg*rhov*mul, 1.e-6), 0.);
    cc = max(7.0e-4 * laml * pow(prl,0.35), 0.0);
    tchf = tsat + pow(qchf,0.3)/max( cc * pow(aa,0.2) * pow(bb,0.7), 1.e-6 );

    if( tsf > tchf ){
        iflg = 1;
        qnb = qchf;
    }
    else{
        iflg = 0;
    // ! Kutateladze boiling curve
        qnb = pow(cc,3.33) * pow(max(tsf-tsat, 0.),3.33) 
            * pow(aa,0.67) * pow(bb,2.33);
        qnb = min(qnb, qchf);
    }
    // end if

    // ! consider void factor
    // !  reduce qnb when alp > 0.3 and make it 0 at alp = 0.95
    // !  with a linear function
    if( alp > 0.3){
        fa = max( 1.-1./0.65*(alp-0.3) , 0. );
        qnb = qnb*fa;
    }
    // end if

    return; 
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void htranspar(double tsf,double dm,double epsrad, double vrv,double vrl,
               double p,double alp, double tv,double tl,double tsat,double hfg,
               double rhov,double rhol,double cpv,double cpl,double muv,
               double mul,double lamv,double laml,double sig,double rhovf,
               double cpvf,double lamvf,double muvf,
               int irad,
               int &icreg, 
               double &qflux, double &qboil, double &qcon, double &qrad)
{
    /* 
    tsf,    & ! melt temperature (surface)
    dm,     & ! melt particle dia
    epsrad, & ! melt emissivity
    vrv,vrl,& ! relative velocity of particle and gas/liquid coolant
    alp,    & ! void fraction
    p,      & ! pressure
    tv,     & ! gas temperature
    tl,     & ! water temperature
    tsat,   & ! water saturation temperature (for boling)
    hfg,    & ! water latent heat for evaporation (for boiling)
    rhov,rhol,muv,mul,lamv,laml,cpv,cpl,sig,& ! water/vapor props.
    rhovf,cpvf,lamvf,muvf ! vapor film props. at (tsat+tsf)/2

    integer,intent(in) :: irad ! flag for radiation effect 
        ! (1:do radiation, 0:no radiation)
    real(kind(1d0)),intent(out) :: qflux, qboil, qcon, qrad   ! heat fluxes
    integer,intent(out) :: icreg  ! heat transfer regime
        !  1: c convection
        !  2: n nucleate boiling
        !  3: f film boiling
        !  4: t transition boiling
        !  5: p convection prevails boiling
    real(kind(1d0)) :: qfb, qnb, qcvv, qcvl, nucvv, nucvl,  &
                    tbini, tminfb, tchf, prv, prl, rev, rel, facv
    integer :: ifb, inb
    character(len=1) :: creg 
    */
    double qfb, qnb, qcvv, qcvl, nucvv, nucvl, tbini, tminfb, tchf, prv, prl, rev, rel, facv;
    int ifb, inb;
    char creg;

    /* ! boiling inception temp.
    ! simply tsf>tsat is used because surface temp. is given.
    ! if tm (inside melt) is given, the condition would be
    !  ti=(bm*tm+bl*tl)/(bm+bl)>tsat -> tm-tsat > bl/bm*(tsat-tl) */

    tbini = tsat ;
    ifb = -1; // ! initialize flags
    inb = -1;

    if (tsf < tbini) { // ! ! low temp. no-boiling regime
        qboil = 0; 
        qrad = 0; 
        creg = 'c';
    } 
    else{ // ! now, boiling is considered
        htfboilpar(p, alp, dm, vrl, tl, tsf, epsrad, tsat, hfg,  
                        rhol, cpl, mul, laml, sig, rhovf,cpvf,lamvf,muvf, 
                        qfb, qrad, tminfb,ifb);
        if (ifb == 0){      // ! film boiling
            qboil = qfb;
            creg = 'f';
        }
        else{                      // ! not film boiling
            htnucboil(p,alp,tl,tsf, tsat, hfg, rhov,rhol,cpl,laml,mul,sig, 
                        qnb,tchf,inb); // ! nuc. boiling model
            if(inb == 0){        // ! nuc. boil.
                qboil = qnb;
                creg = 'n';
            }
            else{                       // ! not nuc. boil. it's transition
                // ! transition boiling { linear interplation between
                // ! (tchf,qchf)--(tminfb,qmfb)
                if(tminfb > tchf){
                    qboil = qnb-(qnb-qfb)/(tminfb-tchf)*(tsf-tchf);
                }
                else{                     // ! this should not happen, but...
                    qboil = 0.5*(qfb+qnb);
                }
                // end if
                creg = 't';
            }
            // end if
        }
        // end if
    }
    // end if 

    // !  convection (Ranz & Marshall [RAN52]) 
    prl = cpl*mul/laml;
    prv = cpv*muv/lamv;
    rel = vrl*dm*rhol/mul;
    rev = vrv*dm*rhov/muv;
    nucvl = min(2.+0.6*pow(rel,0.5) * pow(prl,0.333), 100.);
    qcvl = nucvl*laml/max(dm, 1e-3)*(tsf-tl);
    nucvv = min(2.+0.6*pow(rev,0.5) * pow(prv,0.333), 100.);
    qcvv = nucvv*lamv/max(dm, 1e-3)*(tsf-tv);

    // !  composition of liq. and vap. convections
    if (alp < 0.3){
        facv = 0.;
    }
    else if (alp < 0.75){
        facv = ( alp - 0.3 ) / ( 0.75 - 0.3 );
    }
    else{
        facv = 1.;
    }
    // end if

    qcon = (1.-facv)*qcvl + facv*qcvv;     // ! this is for non-boiling
    if(tsf>tbini && inb==0){   // ! in nuc. boiling regime
        if(qcon > qboil){       // ! convenction prevails boiling
            qboil = 0.;          // ! in nuc. boiling, surface is wettable
            creg = 'p';          // !{ liq. conv. is possible.
        }                       // ! take only conv. if it prevails
        else{                   // ! gas conv. with nuc. boiling
            qcon = facv*qcvv;
        }
        // end if
    }
    else if(tsf>tbini){         // ! gas conv. with trans/film boiling
        qcon = facv*qcvv;       // ! conv. only with gas
                                // ! use max side of qboil and qcvv.
        if(qcon > qboil){       // ! (void is considered in film boil,
            qboil = 0.;          // !  they are exclusive)
            creg = 'p';
        }
        else{
            qcon = 0.;
        }
        // end if
    }
    // end if

    // ! construct total heat flux
    qflux = qboil + qcon + 7./8.*qrad;

    // ! set integer vars for mode
    if(creg == 'c'){
        icreg = 1;}
    else if(creg == 'n'){
        icreg = 2;}
    else if(creg == 'f'){
        icreg = 3;}
    else if(creg == 't'){
        icreg = 4;}
    else if(creg == 'p'){
        icreg = 5;}
    else{
        icreg = 0;}  // ! error 
    // end if

    return;
    // return icreg;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
/* int main ()
{
    double a =1,b=1,c=1,d=1;
    int i1=-1; 
    htranspar(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                i1,a,b,c,d);


    cout << i1 << endl; 

    cout << "DONE"<< endl; 

    return 0;
} */
#endif