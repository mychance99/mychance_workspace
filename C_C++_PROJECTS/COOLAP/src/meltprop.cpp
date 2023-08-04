#ifndef MELTPROP_H
#define MELTPROP_H

// meltprop.cpp
using namespace std;
#include <string>
#include <cmath>


// ############################################################################ // 
// STRUCT & Parameters
// ############################################################################ // 

// ############################################################################ // 
double tempm0base = 300, enem0base;
// ############################################################################ // 

// ############################################################################ // 
// ! data type for basic parameters to calculate melt physical properties (used in unified functions)
struct mpropbase
{
    string mname;
    double tliqu,tsoli,tmelt,  // ! liquidus/solidus temp., average of them
      cpliq, cpsol, lheat, rholiq, rhosol,  // ! spec. heat, latent heat, density
      lamliq, lamsol, visc, stens, emmis; // ! thermal conductivity, viscosity, surface tension, emmisivity
};
// ############################################################################ // 

// ############################################################################ // 
// ! data type for "melt properties" for use in the outside
struct mprop
{
    double tm,  // ! melt temperature
      tmel, tsol, tliq,     // ! melting point (ave(tsol,tliq)), solidus and liquidus points 
      rhom,enem,cpm,dedt,   // ! melt density, internal energy, specific heat,dedt,
      mum,lamm,sigm,epsm;     // ! viscosity, thermal cond., surface tension, emmissivity
    // ! cpm = specific heat as physical properties, not involving latend heat
    // ! latent heat is involved in cp in solidus-liquidus span for calculating int. energy
    // ! and its temperature derivative, dedt
};
// ############################################################################ // 

























// FUNCTIONS



// ############################################################################ // 
void mtempene (mpropbase mpara, double tm, double &enem)
// ! calculate melt energy from temperature
// ! latent heat is considered by embedding it in specific heat
// ! for the solidus-liquidus span
// ! specific heat is given by linear function of temp.
{
    double c1, c2, c3;
    c1=0.5*(4*mpara.lheat/(mpara.tliqu-mpara.tsoli)-(mpara.cpsol-mpara.cpliq)) 
        /(mpara.tliqu-mpara.tsoli);
    c2=-0.5*(4*mpara.lheat/(mpara.tliqu-mpara.tsoli)+(mpara.cpsol-mpara.cpliq))
        /(mpara.tliqu-mpara.tsoli);
    c3=mpara.tmelt*(mpara.cpsol-mpara.cpliq)+mpara.lheat;
    
    if (tm < mpara.tsoli)
        {enem = mpara.cpsol*tm;}
    else if (tm < mpara.tmelt)
        {enem = mpara.cpsol*tm+ c1*pow((tm-mpara.tsoli),2);}
    else if (tm < mpara.tliqu)
        {enem = mpara.cpliq*tm+ c2*pow((tm-mpara.tliqu),2)+c3;}
    else
        {enem = mpara.cpliq*tm+c3;}
    // return enem;
}
// ############################################################################ // 


// ############################################################################ // 
void mdenedtemp (mpropbase mpara, double tm, double &dedt)
// ! temperature derivative of internal energy
{
    double c1, c2;

    c1=0.5*(4*mpara.lheat/(mpara.tliqu-mpara.tsoli)-(mpara.cpsol-mpara.cpliq)) 
        /(mpara.tliqu-mpara.tsoli);
    c2=-0.5*(4*mpara.lheat/(mpara.tliqu-mpara.tsoli)+(mpara.cpsol-mpara.cpliq)) 
        /(mpara.tliqu-mpara.tsoli);

    if (tm < mpara.tsoli)
        {dedt = mpara.cpsol;}
    else if (tm < mpara.tmelt)
        {dedt = mpara.cpsol + 2.0*c1*(tm-mpara.tsoli);}
    else if (tm < mpara.tliqu)
        {dedt = mpara.cpliq + 2.0*c2*(tm-mpara.tliqu);}
    else
        {dedt = mpara.cpliq;}
    // return dedt;
}
// ############################################################################ // 
























// ############################################################################ // 
void initmprop(mpropbase &mpara)
// ! initialize melt property parameters (mpara) for later call of two subroutines below
{ 
    if (mpara.mname == "corium") // ! corium (OECD/IPS39, SERENA UO2-ZrO2 80:20wt%)
    { 
        mpara.tliqu  =2850;   
        mpara.tsoli=2830;
        mpara.cpliq  =565;    
        mpara.cpsol=445;   
        mpara.lheat=0.362e6;
        mpara.rholiq =7960;   
        mpara.rhosol=9430;
        mpara.lamliq =2.88;   
        mpara.lamsol=2.88;
        mpara.visc   =4.23e-3;  
        mpara.stens=0.45;  
        mpara.emmis=0.79;
    }
    else if (mpara.mname == "corium2") // ! corium (OECD/IPS39, SERENA UO2-ZrO2 80:20wt%)
    {
        mpara.tliqu  =2670.0;   
        mpara.tsoli=2420.0;
        mpara.cpliq  =565.0;    
        mpara.cpsol=445.0;   
        mpara.lheat=0.362e6;
        mpara.rholiq =7960.0;   
        mpara.rhosol=9430.0;
        mpara.lamliq =2.88;   
        mpara.lamsol=2.88;
        mpara.visc   =4.23e-3;  
        mpara.stens=0.45;  
        mpara.emmis=0.79;
    }
    else
    {
        printf("meltprop.initmprop");
    }
    mpara.tmelt=0.5*(mpara.tliqu+mpara.tsoli); // ! average of soliduis/liquidus tempes.
    
    // printf ("call mtempene\n");
    // call mtempene(mpara,tempm0base, enem0base)
    mtempene(mpara, tempm0base, enem0base);
}
// ############################################################################ // 

// ############################################################################ // 
void getmprop(mpropbase mpara, mprop &mpp)
{
    mpp.tmel = mpara.tmelt;
    mpp.tsol = mpara.tsoli;
    mpp.tliq = mpara.tliqu;

    // ! thermal conductivity
    if( mpp.tm < mpp.tmel ) 
        {mpp.lamm = mpara.lamsol;}
    else
        {mpp.lamm = mpara.lamliq;}
    // end if

    // ! viscosity
    mpp.mum  = mpara.visc;

    // ! surface tension
    mpp.sigm = mpara.stens;

    // ! radiation emmisivity
    mpp.epsm = mpara.emmis;

    // ! density
    if(mpp.tm < mpp.tsol) 
        {mpp.rhom = mpara.rhosol;}
    else if( mpp.tm > mpp.tliq) 
        {mpp.rhom = mpara.rholiq;}
    else // ! linear interpolation
        mpp.rhom = mpara.rhosol + (mpp.tm-mpp.tsol) 
            *(mpara.rholiq-mpara.rhosol)/(mpp.tliq-mpp.tsol);
    // end if

    // ! specific heat (as a physical property, latent heat not included)
    if( mpp.tm < mpp.tmel ) 
        {mpp.cpm = mpara.cpsol;}
    else
        {mpp.cpm = mpara.cpliq;}
    // end if

    // ! internal energy and its temperature derivative
    // ! (latent heat considered)
    // call mtempene(mpara,mpp->tm, mpp->enem)
    // call mdenedtemp(mpara,mpp->tm, mpp->dedt)
    mtempene(mpara, mpp.tm, mpp.enem);
    mdenedtemp(mpara, mpp.tm, mpp.dedt);
}
// ############################################################################ // 

// ############################################################################ // 
void menetemp (mpropbase mpara, double enem, double &tm)
// ! reverse calculate melt temperature from the energy
{
    double to,tn,eo,dedt,delta;

    // ! initial guess
    to = mpara.tmelt;

    // ! newton conversion
    for (int itr = 1; itr <= 20; itr++) // do itr=1, 20
    {
        
        // call mtempene(mpara,to, eo)
        // call mdenedtemp(mpara,to, dedt)
        mtempene(mpara, to, eo);
        mdenedtemp(mpara, to, dedt);

        tn = to-(eo-enem)/dedt;
        delta = fabs(eo/enem-1);
        if(delta < 1e-5) 
        {
            break;
        }
        to = tn;
    }
    if(delta > 1e-5)
    {
        printf("meltprop.menetemp\n");
    }
        // write(*,fmt='(a,1p,2(1x,e11.3e3))') &
        //   'menetemp: itr>20 no conv.// ! eo,to=', eo,to
    //   end if

    tm = tn;

    // return tm;
}
// ############################################################################ // 

// ############################################################################ // 
// ############################################################################ // 

// ! calculate melt physical properties from given temperature


// ############################################################################ // 
// TMP

// int main ()
// {
//     mpropbase mpara;
//     mpara.mname = "corium2";
//     initmprop(&mpara);

//     mprop mpp;
//     getmprop(mpara, &mpp);

//     mpp.tm = 2000;
//     mpp.enem = mtempene(mpara, mpp.tm, mpp.enem);


//     printf("DONEPRE\n");
//     printf("DONE");
//     return 0;
// }

#endif