#ifndef IDGAS_H
#define IDGAS_H

using namespace std;

#include "consts.cpp"
#include "wrsttab.cpp"
/*###############################################################################*/ 


// struct drpt
// {
//     double drdp, drdt;
// }
/*###############################################################################*/ 
void stigas1(double p,double t,double mm,double cc, 
             double &rho,double &e,double &cv,double &cp)
{
    double t0 = 273.15;

    if(t<=0.0 || p<0.0 || mm<=0.0 || cc<=0.0) { // ! zero p is allowed.
        printf("stigas: Input out of range# ! (range: p>=0; t,mm,cc>0)\n");
        printf("  p(Pa)= %f \n", p);
        printf("  t(K) = %f \n", t);
        printf("  cc   = %f \n", cc);
        printf("  mm   = %f \n", mm);
        exit(EXIT_FAILURE);
    }

    rho = p*mm/rgasgen/t;  // ! rgasgen: ideal gas const. R=8.31(J/molK) in "consts" module
    cv = (cc-1.0)*rgasgen/mm;
    cp = cv + rgasgen/mm;
    e = cv*(t-t0);
}

void stigas2(double p,double t,double mm,double cc, 
             double &rho,double &e,double &cv,double &cp, double &drdp, double &drdt)
{
    double t0 = 273.15;

    if(t<=0.0 || p<0.0 || mm<=0.0 || cc<=0.0) { // ! zero p is allowed.
        printf("stigas: Input out of range# ! (range: p>=0; t,mm,cc>0)\n");
        printf("  p(Pa)= %f \n", p);
        printf("  t(K) = %f \n", t);
        printf("  cc   = %f \n", cc);
        printf("  mm   = %f \n", mm);
        exit(EXIT_FAILURE);
    }

    rho = p*mm/rgasgen/t;  // ! rgasgen: ideal gas const. R=8.31(J/molK) in "consts" module
    cv = (cc-1.0)*rgasgen/mm;
    cp = cv + rgasgen/mm;
    e = cv*(t-t0);

    drdp = mm/rgasgen/t;
    drdt = -p*mm/rgasgen/pow(t,2);
    return;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
/* 
void mixpropgas (t,n,name,p,  mumix,lammix,cpmix,betamix)
// ! evaluate gas mixture viscosity and heat conductivity
// ! if h2o is specified, saturation pressure at given T is assumed.
// ! density is out of the scope, it needs pressures being given.
{

} 
*/
/*###############################################################################*/ 



#endif