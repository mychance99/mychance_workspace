#ifndef INTTAB_H
#define INTTAB_H

#include <iostream>
#include <math.h>
#include <random>
#include <vector>

using namespace std;

/*###############################################################################*/
// STRUCT
const int max_tbl = 500;

struct inttab
{
    int n=0, typ=0;
    // double u[max_tbl]={}, x[max_tbl]={};
    vector<double> u, x;
};

struct timetab
{
    int n=0, itu=0, typ=0;
    // double t[max_tbl]={}, x[max_tbl]={};
    vector<double> t, x;
};
/*###############################################################################*/


/*###############################################################################*/
void itabcheck(inttab tab) 
{
    int i,j;
    if (tab.n < 2){return;}

    for (i=1; i<tab.n; i++) 
    {
        if (tab.u[i-1] > tab.u[i])
        {
            printf("Reversed position specification found in an itab data! \n");
            printf("TMP \n");
            exit(EXIT_FAILURE);
        }
    }
}
/*###############################################################################*/

/*###############################################################################*/
void itabeval(inttab itab,double poi,double &val)
{
    int i;
    if(itab.n == 0){val = itab.x[0];}
    else if(poi <= itab.u[0]){val = itab.x[0];}
    else if(poi >= itab.u[itab.n-1]){val = itab.x[itab.n-1];}
    else {
        for(i=1; i<itab.n; i++)
        {
            if( poi >= itab.u[i-1] && poi <= itab.u[i] )
            { 
                if(itab.u[i-1]-itab.u[i] == 0.0)
                {
                    val = 0.5*(itab.x[i-1] + itab.x[i]);
                }
                else
                {
                    if(itab.typ == 0){  //  ! typ=0 means linear interpolation
                        val=itab.x[i-1] 
                            +(itab.x[i]-itab.x[i-1])/(itab.u[i]-itab.u[i-1]) 
                            *(poi-itab.u[i-1]);
                    }
                    else{                    //  ! typ=1 means stepwise evolution
                        val=itab.x[i-1];
                    }
                }
            }
        }
    }
}
/*###############################################################################*/

/*###############################################################################*/
void itabevalr (inttab itab,double val,double &poi)
{
    int i;
    if(itab.n == 0){poi = itab.u[0];}
    else if(val <= itab.x[0]){poi = itab.u[0];}
    else if(val >= itab.x[itab.n-1]){poi = itab.u[itab.n-1];}
    else {
        for(i=1; i<itab.n; i++)
        {
            if( val >= itab.x[i-1] && val <= itab.x[i] )
            { 
                if(itab.x[i-1]-itab.x[i] == 0.0)
                {
                    poi = 0.5*(itab.u[i-1] + itab.u[i]);
                }
                else
                {
                    if(itab.typ == 0){  //  ! typ=0 means linear interpolation
                        poi=itab.u[i-1] 
                            +(itab.u[i]-itab.u[i-1])/(itab.x[i]-itab.x[i-1]) 
                            *(val-itab.x[i-1]);
                    }
                    else{                    //  ! typ=1 means stepwise evolution
                        poi=0.5*(itab.u[i-1]+itab.u[i]);
                    }
                }
            }
        }
    }
}
/*###############################################################################*/

/*###############################################################################*/
void ttabcheck(timetab tab)
{
    int i, j;
    if(tab.n <2){return;}
    for(i=1; i<tab.n; i++)
    {
        if(tab.t[i-1] > tab.t[i])
        {
            printf ("Reversed time specification found in a ttab data! \n");
            printf ("tmp\n");
            // # print ("Time : ", tab.t[j], j=i, i+1)
            exit(EXIT_FAILURE);
        }
    }
}
/*###############################################################################*/

/*###############################################################################*/
void ttabconv(timetab ttab)
{
    printf("NEED TO CHECK - ttabconv \n");
    // if(ttab.itu == 1){ttab.t[0:ttab.n-1] = ttab.t[0:ttab.n-1]*3600.0;}
    // else if(ttab.itu == 2){ttab.t[0:ttab.n-1] = ttab.t[0:ttab.n-1]*3600.0*24.0;}
    ttab.itu = 0;
}
/*###############################################################################*/

/*###############################################################################*/
void ttabeval(timetab ttab, double time, double &val)
{
    int i;
    if(ttab.itu != 0){ttabconv(ttab);}

    if(ttab.n == 0){val = ttab.x[0];}
    else if(time <= ttab.t[0]){val = ttab.x[0];}
    else if(time >= ttab.t[ttab.n-1]){val = ttab.x[ttab.n-1];}
    else {
        for(i=1; i<ttab.n; i++)
        {
            if( time >= ttab.t[i-1] && time <= ttab.t[i] )
            { 
                if(ttab.t[i-1]-ttab.t[i] == 0.0)
                {
                    val = 0.5*(ttab.x[i-1] + ttab.x[i]);
                }
                else
                {
                    if(ttab.typ == 0){  //  ! typ=0 means linear interpolation
                        val=ttab.x[i-1] 
                            +(ttab.x[i]-ttab.x[i-1])/(ttab.t[i]-ttab.t[i-1]) 
                            *(time-ttab.t[i-1]);
                    }
                    else{                    //  ! typ=1 means stepwise evolution
                        val=ttab.x[i-1];
                    }
                }
            }
        }
    }
}
/*###############################################################################*/

/*###############################################################################*/
void showitab(inttab tab, int unt=0)
{
    int i, u=0;
    printf("LATER -showitab \n");
}
/*###############################################################################*/

/*###############################################################################*/
void showttab(timetab tab, int unt=0)
{
    printf("LATER -showttab \n");
}
/*###############################################################################*/

/*###############################################################################*/
// int main()
// {
//     struct inttab tab; 
//     tab.x[0]=44;
//     printf("hello \n");

//     // // melthex mh;
//     // // mh.nthex = 3;
//     // meltdebr md;
//     // float time=1, dt=0.1;
//     // // std::cout << typeid(dt).name() << std::endl;
//     // md.psieve.nn=20;
//     // printf ("\n %d \n", md.psieve.nn);
//     // // mjetevol(&md, time, dt);
//     // mjetevol(md, time, dt);
//     // printf ("\n %d \n", md.psieve.nn);
    

//     printf(" ################\n ");
//     printf("DONE\n");

//     return 0;
// }
#endif