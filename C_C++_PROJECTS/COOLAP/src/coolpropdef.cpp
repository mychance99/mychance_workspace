#ifndef COOLPROPDEF_H
#define COOLPROPDEF_H

using namespace std;
#include <string>


const int nmax = 10;

// ! data type for "coolant properties" as ambient fluids of heat transfer body
// ! assuming gas mixture and water in thermal equilibrium (at steam partial pressure)
// ! some thermodynamic properties (enthalpy or internal energy, hfg) are not included

struct coolprop
{
    int ng=0;
    string gname[nmax];
    double ptot,pg[nmax], tl,tg,tspt, // ! pressures (total,steam  noncondensibles)
                                      // !  temps. (tspt = saturation temp. at ptot)
          hfg,  // ! latent heat of water boiling (hfg(tspt))
          rhol,rhov,rhog,cpl,cpv,cpg,mul,muv,mug,laml,lamv,lamg, sig,
          betal,betav,betag,         // ! water(l), steam(v) and gas mixture(g) props.
                                      // ! steam properties are at ptot (for boiling heat transfer)
          hlev, wvol, wmdi;            // ! water level, volume, flow-in rate in the can
};



// int main ()
// {
//     coolprop clp; 
//     printf("DNE\n");
//     return 0;
// }

#endif
