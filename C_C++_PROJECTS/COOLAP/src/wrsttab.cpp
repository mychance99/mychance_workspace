#ifndef WRSTTAB_H
#define WRSTTAB_H

#include <iostream>
#include <math.h>
using namespace std;

/*###############################################################################*/ 
const int mp1=11, mt1g=12+1, mt1l=10+1; // 10 12 10
int jmaxg[mp1]={}, jmaxl[mp1]={};
double p1mesh[mp1],t1meshg[mt1g],t1meshl[mt1l], 
       tabdeng[mp1][mt1g][3], tabdenl[mp1][mt1l][3], 
       tabentg[mp1][mt1g][3], tabentl[mp1][mt1l][3];
const double pcwater = 22.12e6, tcwater = 647.3;
/*###############################################################################*/ 

/*###############################################################################*/ 
void ipcub1d(double x0,double f0,double fd0, double x1,double f1,double fd1,double x,
             double &f,double &fd,int &i)
{
    /* ! cubic interpolation 1d */
    // ! i : flag (0:good, -1:error)
    // ! x0,x1 : edges of the part
    // ! f0,fd0, f1,fd1 : f and f' at x0,x1
    // ! x : place to interpolate
    // ! f,fd : f and f' at x (output)
    double  a,b, phi,del;

    del = x1-x0;
    if(del == 0.0)
    {
        f = 0.0;
        fd = 0.0;
        i = -1;
        return;
    }
    // endif;
    // !;
    phi = (f1-f0)/del;
    // a = (-2.0*phi + fd1)/del/del  + fd0/del/del;
    a = (-2.0*phi + fd1 + fd0)/del/del;
    b = ( 3.0*phi - fd1 - 2.0*fd0)/del;
    // !;
    f = a*pow((x-x0),3) + b*pow((x-x0),2) + fd0*(x-x0) +f0;
    fd = 3.0*a*pow((x-x0),2) + 2.0*b*(x-x0) + fd0;
    i = 0;
    return;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void iplin1d(double x0,double f0, double x1,double f1, double x, 
             double &f,double &fd,int &i)
{
    /* ! linear interpolation 1d */
    // ! f,fd : f and f' at x (output)
    // ! i : flag (0:good, -1:error)
    // ! x0,x1 : edges of the part
    // ! f0, f1 : f at x0,x1
    // ! x : place to interpolate
    // ! f : f at x (output)
    // ! fd : gradient of f
    double a, del; 

    del = x1-x0;

    if(del == 0.0)
    {
        f = 0.0;
        fd = 0.0;
        i = -1;
        return;
    }
    a = (f1-f0)/del;
    f = a*(x-x0) + f0;
    fd = a;
    i = 0;

    return;
}
/*###############################################################################*/ 

/* WHANG iptab & eptab -> iptabl, iptabg & eptabl, eptabg */
/*###############################################################################*/ 
void iptabl(double x,double y,int mx,int my,
            int jmax[],double tab[][mt1l][3],double xmesh[],double ymesh[], 
            double &f,double &fx,double &fy,int &kx,int &ky,int &ifx,int &ify)
{
    // ! local vars
    int i,j, itmp, i1;
    double  x0,x1,y0,y1, f00,f10,f01,f11, 
         fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11, 
         f0a,fx0a,fy0a, f1a,fx1a,fy1a, c0a, c1a, del;

    // ! search for node including point (x,y)
    if (x < xmesh[0])
    {
        ifx = -1;
        i = 0;
    }
    else if(x > xmesh[mx-1])
    {
        ifx = +1;
        i = mx;
    }
    else
    {
    // !mori021028 make robust even if x is NaN
        if(x >= xmesh[0] && x <= xmesh[mx-1])
        {
            ifx = 0;
            i=0;
            // 100       continue
            i=i+1;
            while (x > xmesh[i]) {i+=1;}
            // if(x .gt. xmesh(i)) goto 100
        }
        else
        {
            ifx = -9;
            i=0;
        }
        // end if
    }
    // end if
    // ! now, i shows the box in which x is.
    // ! ymesh direction depends on gas/liq

    // ! liq
    if(ymesh[my-1] > ymesh[0])
    {
        if(ifx == 1){i1 = mx+1;}
        else if(ifx == -1){i1 = 1;}
        else {i1 = i;}
        // end if
        // !
        if(y < ymesh[0])
        {
            ify = -1;
            j = 0;
        }
        else if(y > ymesh[jmax[i1-1]])
        {
            ify = +1;
            j = jmax[i1-1];
        }
        else
        // !mori021028 make robust even if y is NaN
        {
            if(y >= ymesh[0] && y <= ymesh[jmax[i1-1]])
            {
                ify = 0;
                j=0;
                // 200          continue
                j=j+1;
                while (y > ymesh[j]){j+=1;}
                //  if(y .gt. ymesh(j)) goto 200
            }
            else
            {
                ify = -9;
                j=0;
            }
            //  end if
        }
        //  end if
    }
    // ! gas
    else 
    {
        i1 = i;
    // !
        if(y < ymesh[jmax[i1]])
        {
            ify = -1;
            j = jmax[i1];
        }
        else if(y > ymesh[0])
        {
            ify = +1;
            j = 0;
        }
        else
        // !mori021028 make robust even if y is NaN
        {
            if(y >= ymesh[jmax[i1]] && y <= ymesh[0])
            {
                ify = 0;
                j=0;
            // 210          continue
                j=j+1;
                while (y < ymesh[j]) {j += 1;}
            // if(y .lt. ymesh(j)) goto 210
            }
            else
            {
                ify = -9;
                j=0;
            }
            // end if
        }
        // end if
    }
    // end if

    kx = i;
    ky = j;

    // ! out-of-range handling;
    if((ifx != 0) || (ify != 0))
    {
        f = 0.0;
        fx = 0.0;
        fy = 0.0;
        return;
    }
    // endif;

    // ! interpolation;
    x0 = xmesh[i-1];
    x1 = xmesh[i  ];
    y0 = ymesh[j-1];
    y1 = ymesh[j  ];
    f00 = tab[i-1][j-1][0];
    f10 = tab[i  ][j-1][0];
    f01 = tab[i-1][j  ][0];
    f11 = tab[i  ][j  ][0];
    fx00 = tab[i-1][j-1][1];
    fx10 = tab[i  ][j-1][1];
    fx01 = tab[i-1][j  ][1];
    fx11 = tab[i  ][j  ][1];
    fy00 = tab[i-1][j-1][2];
    fy10 = tab[i  ][j-1][2];
    fy01 = tab[i-1][j  ][2];
    fy11 = tab[i  ][j  ][2];

    // ! vertical -> horizontal interpolation
    // !
    // !  01               11     1) cubic  int. on 00-01 line for f0a, fy0a
    // !    +-------------+       2) linear int. on 00-01 line for fx0a, c0a 
    // !    |             |           (c0a: gradient of fx on 00-01)
    // !    |             |       3) do same as 1 and 2 on 10-11 line for
    // ! 0a *--------X----* 1a       f1a, fy1a, fx1a, c1a
    // !    |             |       4) cubic int. on 0a-1a line for f and fx
    // !    |             |       5) get fy by special interpolation method
    // !    +-------------+          derived considering variance of f(x)
    // !  00               10        for small change of y

    // ! on 00-01 line
    ipcub1d(y0,f00,fy00, y1,f01,fy01, y, f0a, fy0a, itmp);
    iplin1d(y0,fx00,y1,fx01, y, fx0a, c0a, itmp);
    // ! on 10-11 line
    ipcub1d(y0,f10,fy10, y1,f11,fy11, y, f1a,fy1a, itmp);
    iplin1d(y0,fx10,y1,fx11,y, fx1a, c1a, itmp);
    // ! on 0a-1a line
    ipcub1d(x0,f0a,fx0a, x1,f1a,fx1a, x, f,fx, itmp);

    del = x1-x0;

    fy = (-2.0/pow(del,3)*(fy1a-fy0a)+1.0/pow(del,2)*(c0a+c1a))*pow((x-x0),3)
        +( 3.0/pow(del,2)*(fy1a-fy0a)-1.0/del*(2.0*c0a+c1a))*pow((x-x0),2)
        +c0a*(x-x0) + fy0a;

    return; 
}
void iptabg(double x,double y,int mx,int my,
            int jmax[],double tab[][mt1g][3],double xmesh[],double ymesh[], 
            double &f,double &fx,double &fy,int &kx,int &ky,int &ifx,int &ify)
{
    // ! local vars
    int i,j, itmp, i1;
    double  x0,x1,y0,y1, f00,f10,f01,f11, 
         fx00,fx10,fx01,fx11,fy00,fy10,fy01,fy11, 
         f0a,fx0a,fy0a, f1a,fx1a,fy1a, c0a, c1a, del;

    // ! search for node including point (x,y)
    if (x < xmesh[0])
    {
        ifx = -1;
        i = 0;
    }
    else if(x > xmesh[mx-1])
    {
        ifx = +1;
        i = mx;
    }
    else
    {
    // !mori021028 make robust even if x is NaN
        if(x >= xmesh[0] && x <= xmesh[mx-1])
        {
            ifx = 0;
            i=0;
            // 100       continue
            i=i+1;
            while (x > xmesh[i]) {i+=1;}
            // if(x .gt. xmesh(i)) goto 100
        }
        else
        {
            ifx = -9;
            i=0;
        }
        // end if
    }
    // end if
    // ! now, i shows the box in which x is.
    // ! ymesh direction depends on gas/liq

    // ! liq
    if(ymesh[my-1] > ymesh[0])
    {
        if(ifx == +1){i1 = mx+1;}
        else if(ifx == -1){i1 = 1;}
        else {i1 = i;}
        // end if
        // !
        if(y < ymesh[0])
        {
            ify = -1;
            j = 0;
        }
        else if(y > ymesh[jmax[i1-1]])
        {
            ify = +1;
            j = jmax[i1-1];
        }
        else
        // !mori021028 make robust even if y is NaN
        {
            if(y >= ymesh[0] && y <= ymesh[jmax[i1-1]])
            {
                ify = 0;
                j=0;
                // 200          continue
                j=j+1;
                while (y > ymesh[j]){j+=1;}
                //  if(y .gt. ymesh(j)) goto 200
            }
            else
            {
                ify = -9;
                j=0;
            }
            //  end if
        }
        //  end if
    }
    // ! gas
    else 
    {
        i1 = i;
    // !
        if(y < ymesh[jmax[i1]])
        {
            ify = -1;
            j = jmax[i1];
        }
        else if(y > ymesh[0])
        {
            ify = +1;
            j = 0;
        }
        else
        // !mori021028 make robust even if y is NaN
        {
            if(y >= ymesh[jmax[i1]] && y <= ymesh[0])
            {
                ify = 0;
                j=0;
            // 210          continue
                j=j+1;
                while (y < ymesh[j]) {j += 1;}
            // if(y .lt. ymesh(j)) goto 210
            }
            else
            {
                ify = -9;
                j=0;
            }
            // end if
        }
        // end if
    }
    // end if

    kx = i;
    ky = j;

    // ! out-of-range handling;
    if((ifx != 0) || (ify != 0))
    {
        f = 0.0;
        fx = 0.0;
        fy = 0.0;
        return;
    }
    // endif;

    // ! interpolation;
    x0 = xmesh[i-1];
    x1 = xmesh[i  ];
    y0 = ymesh[j-1];
    y1 = ymesh[j  ];
    f00 = tab[i-1][j-1][0];
    f10 = tab[i  ][j-1][0];
    f01 = tab[i-1][j  ][0];
    f11 = tab[i  ][j  ][0];
    fx00 = tab[i-1][j-1][1];
    fx10 = tab[i  ][j-1][1];
    fx01 = tab[i-1][j  ][1];
    fx11 = tab[i  ][j  ][1];
    fy00 = tab[i-1][j-1][2];
    fy10 = tab[i  ][j-1][2];
    fy01 = tab[i-1][j  ][2];
    fy11 = tab[i  ][j  ][2];

    // ! vertical -> horizontal interpolation
    // !
    // !  01               11     1) cubic  int. on 00-01 line for f0a, fy0a
    // !    +-------------+       2) linear int. on 00-01 line for fx0a, c0a 
    // !    |             |           (c0a: gradient of fx on 00-01)
    // !    |             |       3) do same as 1 and 2 on 10-11 line for
    // ! 0a *--------X----* 1a       f1a, fy1a, fx1a, c1a
    // !    |             |       4) cubic int. on 0a-1a line for f and fx
    // !    |             |       5) get fy by special interpolation method
    // !    +-------------+          derived considering variance of f(x)
    // !  00               10        for small change of y

    // ! on 00-01 line
    ipcub1d(y0,f00,fy00, y1,f01,fy01, y, f0a, fy0a, itmp);
    iplin1d(y0,fx00,y1,fx01, y, fx0a, c0a, itmp);
    // ! on 10-11 line
    ipcub1d(y0,f10,fy10, y1,f11,fy11, y, f1a,fy1a, itmp);
    iplin1d(y0,fx10,y1,fx11,y, fx1a, c1a, itmp);
    // ! on 0a-1a line
    ipcub1d(x0,f0a,fx0a, x1,f1a,fx1a, x, f,fx, itmp);

    del = x1-x0;

    fy = (-2.0/pow(del,3)*(fy1a-fy0a)+1.0/pow(del,2)*(c0a+c1a))*pow((x-x0),3)
        +( 3.0/pow(del,2)*(fy1a-fy0a)-1.0/del*(2.0*c0a+c1a))*pow((x-x0),2)
        +c0a*(x-x0) + fy0a;

    return; 
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void epidgr(double p10,double t10, double f0,double fp10,double ft10,double p1,
            double t1, int imode, 
            double &f,double &fp1,double &ft1)
{
    double p0,t0,fp0,ft0, p,t,fp,ft, pc,tc, a;

    pc = 22.12;
    tc = 647.3;

    // ! p10, t10, f0, fp10, ft10 input (the base p1, t1, f and derivatives)
    // ! p1, t1           input (p1, t1 for extrapolation) 
    // ! imode            0:t>tmax, 1:p<pmin, 2:t>tmax and p<pmin
    // ! f, fp1, ft1      output (extrapolation result)
    // !  (i/o variables are all normalized (log10(p/pc), t/tc-1))

    p0 = pow(10.0,p10) *pc;
    p  = pow(10.0,p1)  *pc;
    t0 = (t10+1.0) *tc;
    t  = (t1 +1.0) *tc;
    fp0 = fp10 /p0/log(10.0);
    ft0 = ft10 /tc;

    // ! simple ideal gas law (T derivative not smooth)
    // !  f=a(p/t)
    a = f0*t0/p0;
    f = a*p/t;
    // !
    // ! derivative differs for extrapolation modes
    // !      
    if(imode == 0)
    {
        fp = fp0*t0/t;
        ft = -a*p/pow(t,2);
    }
    else if(imode == 1)
    {
        fp = a/t;
        ft = ft0*p/p0;
    }
    else
    {
        fp = a/t;
        ft = -a*p/pow(t,2);
    }
    //end if
    // !
    fp1 = fp*p*log(10.0);
    ft1 = ft*tc;

    return;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void epidgh(double p10,double t10, double f0,double fp10,double ft10, double ft1p10,
            double p1, double t1, int imode, 
            double &f,double &fp1,double &ft1)
{
    /* 
    ! ideal gas enthalpy law for extrapolation
    
    ! p10, t10, f0, fp10, ft10 input (the base p1, t1, f and derivatives)
    ! ft1p10 = @ft10/@p1   input
    ! p1, t1           input (p1, t1 for extrapolation) 
    ! imode            0:t>tmax, 1:p<pmin, 2:t>tmax and p<pmin
    ! f, fp1, ft1      output (extrapolation result)
    !  (i/o variables are all normalized (log10(p/pc), t/tc-1)) 
    */
    double p0,t0,fp0,ft0, p,t,fp,ft, pc,tc, a, ftp0;
    pc = 22.12;
    tc = 647.3;

    p0 = pow(10.0,p10) *pc;
    p  = pow(10.0,p1)  *pc;
    t0 = (t10+1.0) *tc;
    t  = (t1 +1.0) *tc;
    fp0 = fp10 /p0/log(10.0);
    ft0 = ft10 /tc;
    ftp0 = ft1p10/tc/p0/log(10.0);

    // !
    // ! simple ideal gas law (T derivative not smooth)
    // !  f= a*(t-t0) + f0
    // !
    a = ft0;
    f = a*(t-t0)+f0;
    // !
    // ! derivative differs for extrapolation modes
    // !      
    if(imode == 0)
    {fp = fp0 + ftp0*(t-t0);}
    else
    {fp = 0.0;}
    // end if
    ft = ft0;
    // !
    fp1 = fp*p*log(10.0);
    ft1 = ft*tc;

    return;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void eplinear(double p10,double t10, double f0,double fp10,double ft10,double fpt10,
            double p1,double t1, int imode, 
            double &f,double &fp1,double &ft1)
{
    /* 
    ! linear law for extrapolation

    ! p10, t10, f0, fp10, ft10 input (the base p1, t1, f and derivatives)
    ! fpt10 = @fp10/@t1 or @ft10/@p1 input
    ! p1, t1           input (p1, t1 for extrapolation) 
    ! imode            0:extrapolate on t (p in the range)
    !                  1:extrapolate on p (t in the range)
    !                  2:extrapolate on both p and t from one point (p0,t0)
    ! f, fp1, ft1      output (extrapolation result)
    !  (i/o variables are all normalized (log10(p/pc), t/tc-1))
    !
    !
    ! f assumed to be proportional to p and t
    ! f  = fp0*(p-p0)+ft0*(t-t0)+f0 
    */
   double p0,t0,fp0,ft0, p,t,fp,ft, pc,tc, a,b, fpt0, ftp0;

    pc = 22.12;
    tc = 647.3;

    p0 = pow(10.0,p10) *pc;
    p  = pow(10.0,p1)  *pc;
    t0 = (t10+1.0) *tc;
    t  = (t1 +1.0) *tc;
    fp0 = fp10 /p0/log(10.0);
    ft0 = ft10 /tc;
    fpt0 = fpt10/p0/log(10.0)/tc;
    ftp0 = fpt0;

    a = fp0;
    b = ft0;
    f = a*(p-p0)+b*(t-t0)+f0;
    fp = a;
    ft = b;
    
    if(imode == 0)
    {fp = fp + ftp0*(t-t0);}
    else if(imode == 1)
    {ft = ft + fpt0*(p-p0);}
    // end if;
    fp1 = fp*p*log(10.0);
    ft1 = ft*tc;

    return; 
}
/*###############################################################################*/ 
// WHANG iptab & eptab -> iptabl, iptabg & eptabl, eptabg
/*###############################################################################*/ 
void eptabl(double x,double y,int kx,int ky,int ifx,int ify,int mx,int my,
           int jmax[],double tab[][mt1l][3],double xmesh[],double ymesh[], int irhg, 
           double &f,double &fx,double &fy)
{
    /* 
    ! extrapolate table (for the input out of table range)

    ! kx,ky: box index returned by iptab()
    ! ifx, ify: output indicating out-of-range status from iptab()
    ! irhg: set 1 for density of gas, 2 for enthalpy of gas, otherwise set 0
    ! x,y: p1, t1
    ! mx,my: max mesh index for p1 and t1
    ! jmax : max valid mesh in t1-direction
    ! tab : table data array
    ! xmesh,ymesh: p1mesh, t1mesh
    ! f,fx,fy: extrapolated result (output) 
    */
    int itmp1; 
    double p10,t10,f0,fp10,ft10,fp1t10,ft1p10;
    // ! p>pmax 
    // !    water: linear
    // !    steam: extrapolate with the ideal gas law

    if (ifx == 1)
    {
        if (ify != 0)
        {
            p10 = xmesh[kx];
            t10 = ymesh[ky];
            f0 = tab[kx][ky][0];
            fp10 = tab[kx][ky][1];
            ft10 = tab[kx][ky][2];
            // ! ep. on pmax edge => get vars. at p10,t1 => then go to p1,t1
            if(irhg == 1)         // ! gas rho
            {epidgr(p10,t10,f0,fp10,ft10, x,y, 2, f, fx, fy);}
            else if(irhg == 2)    // ! gas h
            {epidgh(p10,t10,f0,fp10,ft10, 0.0, x,y, 2, f, fx, fy);}
            else                         // ! other (linear)
            {eplinear(p10,t10, f0,fp10,ft10, 0.0, x,y,2, f,fx,fy);}
            // end if
        }
        else
        {
            // ! ip. on pmax edge => get vars. at p10,t1 => then go to p1,t1
            p10 = xmesh[kx];
            t10 = y;
            ipcub1d(ymesh[ky-1],tab[kx][ky-1][0],tab[kx][ky-1][2], 
                                            ymesh[ky],tab[kx][ky][0],tab[kx][ky][2],t10, f0 ,ft10, itmp1);
            iplin1d(ymesh[ky-1],tab[kx][ky-1][1], ymesh[ky],tab[kx][ky][1],
                            t10,fp10,fp1t10, itmp1);

            if(irhg == 1)
            {epidgr(p10,t10, f0,fp10,ft10, x,y, 1,  f,fx,fy);}
            else if(irhg == 2)
            {epidgh(p10,t10, f0,fp10,ft10,0.0, x,y, 1, f,fx,fy);}
            else
            {eplinear(p10,t10, f0,fp10,ft10, fp1t10, x,y,1,  f,fx,fy);}
            // end if

        }
    }

    // ! t>tmax or t<tmin, p in the range
    else if (ify != 0 && ifx == 0){
        p10 = x;
        t10 = ymesh[ky];
        ipcub1d(xmesh[kx-1],tab[kx-1][ky][0],tab[kx-1][ky][1], 
            xmesh[kx],tab[kx][ky][0],tab[kx][ky][1],
            p10,f0,fp10,itmp1);
        iplin1d(xmesh[kx-1],tab[kx-1][ky][2],
                xmesh[kx],tab[kx][ky][2], 
                p10,ft10,ft1p10, itmp1);
        if(irhg == 1 && ify == 1) {
            epidgr(p10,t10, f0,fp10,ft10, x,y, 0,  f,fx,fy);
        }
        else if(irhg == 2 && ify == 1) {
            epidgh(p10,t10, f0,fp10,ft10,ft1p10, x,y, 0, 
                f,fx,fy);
        }
        else{
             eplinear(p10,t10, f0,fp10,ft10, ft1p10, 
                x,y,0,  f,fx,fy);
        }
    }
    // ! t in the range, p<pmin
    else if (ify == 0 && ifx == -1) {
       p10 = xmesh[kx];
       t10 = y;
       ipcub1d(ymesh[ky-1],tab[kx][ky-1][0],tab[kx][ky-1][2], 
            ymesh[ky],tab[kx][ky][0],tab[kx][ky][2], 
            t10,f0,ft10,itmp1);
       iplin1d(ymesh[ky-1],tab[kx][ky-1][1], 
            ymesh[ky],tab[kx][ky][1], 
            t10,fp10,fp1t10, itmp1);
       if(irhg == 1) { // ! gas rho
            epidgr(p10,t10,f0,fp10,ft10, x,y, 1, f, fx, fy);
       }   
       else if(irhg == 2) { // ! gas h
            epidgh(p10,t10,f0,fp10,ft10, 0.0, x,y, 1, f, fx, fy);
       } 
       else { // ! other (linear)
            eplinear(p10,t10, f0,fp10,ft10, fp1t10, 
               x,y,1,   f, fx, fy);
       }
    }
    // ! (t>tmax and p<pmin) or (t<tmin and p<pmin) 
    else if (ify != 0 && ifx == -1) {
       p10 = xmesh[kx];
       t10 = ymesh[ky];
       f0 = tab[kx][ky][0];
       fp10 = tab[kx][ky][1];
       ft10 = tab[kx][ky][2];
       if(irhg == 1 && ify == 1){   // ! gas rho
          epidgr(p10,t10,f0,fp10,ft10, x,y, 2, f, fx, fy);
       }
       else if(irhg == 2 && ify == 1){ // ! gas h
          epidgh(p10,t10,f0,fp10,ft10,0.0, x,y, 2, f, fx, fy);
       }
       else {                  // ! other (linear)
          eplinear(p10,t10, f0,fp10,ft10, 0.0, 
               x,y,2,  f, fx, fy);
       }
    }
    //
    else{
        printf("eptab ERROR1\n");
        exit(EXIT_FAILURE);
    }
    printf("eptab.test\n");

    return;
}
void eptabg(double x,double y,int kx,int ky,int ifx,int ify,int mx,int my,
           int jmax[],double tab[][mt1g][3],double xmesh[],double ymesh[], int irhg, 
           double &f,double &fx,double &fy)
{
    /* 
    ! extrapolate table (for the input out of table range)

    ! kx,ky: box index returned by iptab()
    ! ifx, ify: output indicating out-of-range status from iptab()
    ! irhg: set 1 for density of gas, 2 for enthalpy of gas, otherwise set 0
    ! x,y: p1, t1
    ! mx,my: max mesh index for p1 and t1
    ! jmax : max valid mesh in t1-direction
    ! tab : table data array
    ! xmesh,ymesh: p1mesh, t1mesh
    ! f,fx,fy: extrapolated result (output) 
    */
    int itmp1; 
    double p10,t10,f0,fp10,ft10,fp1t10,ft1p10;
    // ! p>pmax 
    // !    water: linear
    // !    steam: extrapolate with the ideal gas law

    if (ifx == 1)
    {
        if (ify != 0)
        {
            p10 = xmesh[kx];
            t10 = ymesh[ky];
            f0 = tab[kx][ky][0];
            fp10 = tab[kx][ky][1];
            ft10 = tab[kx][ky][2];
            // ! ep. on pmax edge => get vars. at p10,t1 => then go to p1,t1
            if(irhg == 1)         // ! gas rho
                {epidgr(p10,t10,f0,fp10,ft10, x,y, 2, f, fx, fy);}
            else if(irhg == 2)    // ! gas h
                {epidgh(p10,t10,f0,fp10,ft10, 0.0, x,y, 2, f, fx, fy);}
            else                         // ! other (linear)
                {eplinear(p10,t10, f0,fp10,ft10, 0.0, x,y,2, f,fx,fy);}
            // end if
        }
        else
        {
            // ! ip. on pmax edge => get vars. at p10,t1 => then go to p1,t1
            p10 = xmesh[kx];
            t10 = y;
            ipcub1d(ymesh[ky-1],tab[kx][ky-1][0],tab[kx][ky-1][2], 
                                            ymesh[ky],tab[kx][ky][0],tab[kx][ky][2],t10, f0 ,ft10, itmp1);
            iplin1d(ymesh[ky-1],tab[kx][ky-1][1], ymesh[ky],tab[kx][ky][1],
                            t10,fp10,fp1t10, itmp1);

            if(irhg == 1)
            {epidgr(p10,t10, f0,fp10,ft10, x,y, 1,  f,fx,fy);}
            else if(irhg == 2)
            {epidgh(p10,t10, f0,fp10,ft10,0.0, x,y, 1, f,fx,fy);}
            else
            {eplinear(p10,t10, f0,fp10,ft10, fp1t10, x,y,1,  f,fx,fy);}
            // end if

        }
    }
    // ! t>tmax or t<tmin, p in the range
    else if (ify != 0 && ifx == 0){
        p10 = x;
        t10 = ymesh[ky];
        ipcub1d(xmesh[kx-1],tab[kx-1][ky][0],tab[kx-1][ky][1], 
            xmesh[kx],tab[kx][ky][0],tab[kx][ky][1],
            p10,f0,fp10,itmp1);
        iplin1d(xmesh[kx-1],tab[kx-1][ky][2],
                xmesh[kx],tab[kx][ky][2], 
                p10,ft10,ft1p10, itmp1);
        if(irhg == 1 && ify == 1) {
            epidgr(p10,t10, f0,fp10,ft10, x,y, 0,  f,fx,fy);
        }
        else if(irhg == 2 && ify == 1) {
            epidgh(p10,t10, f0,fp10,ft10,ft1p10, x,y, 0, 
                f,fx,fy);
        }
        else{
             eplinear(p10,t10, f0,fp10,ft10, ft1p10, 
                x,y,0,  f,fx,fy);
        }
    }
    // ! t in the range, p<pmin
    else if (ify == 0 && ifx == -1) {
       p10 = xmesh[kx];
       t10 = y;
       ipcub1d(ymesh[ky-1],tab[kx][ky-1][0],tab[kx][ky-1][2], 
            ymesh[ky],tab[kx][ky][0],tab[kx][ky][2], 
            t10,f0,ft10,itmp1);
       iplin1d(ymesh[ky-1],tab[kx][ky-1][1], 
            ymesh[ky],tab[kx][ky][1], 
            t10,fp10,fp1t10, itmp1);
       if(irhg == 1) { // ! gas rho
            epidgr(p10,t10,f0,fp10,ft10, x,y, 1, f, fx, fy);
       }   
       else if(irhg == 2) { // ! gas h
            epidgh(p10,t10,f0,fp10,ft10, 0.0, x,y, 1, f, fx, fy);
       } 
       else { // ! other (linear)
            eplinear(p10,t10, f0,fp10,ft10, fp1t10, 
               x,y,1,   f, fx, fy);
       }
    }
    // ! (t>tmax and p<pmin) or (t<tmin and p<pmin) 
    else if (ify != 0 && ifx == -1) {
       p10 = xmesh[kx];
       t10 = ymesh[ky];
       f0 = tab[kx][ky][0];
       fp10 = tab[kx][ky][1];
       ft10 = tab[kx][ky][2];
       if(irhg == 1 && ify == 1){   // ! gas rho
          epidgr(p10,t10,f0,fp10,ft10, x,y, 2, f, fx, fy);
       }
       else if(irhg == 2 && ify == 1){ // ! gas h
          epidgh(p10,t10,f0,fp10,ft10,0.0, x,y, 2, f, fx, fy);
       }
       else {                  // ! other (linear)
          eplinear(p10,t10, f0,fp10,ft10, 0.0, 
               x,y,2,  f, fx, fy);
       }
    }
    //
    else{
        printf("eptab ERROR2\n");
        exit(EXIT_FAILURE);
    }

    return;
}
/*###############################################################################*/ 


/*###############################################################################*/ 
void wrsteamtab (double p,double t,int iphase, 
    double &r,double &drdp,double &drdt, double &h,double &dhdp,double &dhdt,
    int &irng) // irng : optional out
{
    int ifp, ift, kp, kt; 
    double f,fx,fy, p1,t1;
    /* 
    ! iphase: phase index (0:water, 1:steam)
    ! p,t : pressure(Pa), temperature(K)
    ! r, drdp,drdt: density(kg/m3), its derivatives with p and t (output)
    ! h, dhdp,dhdt: enthalpy(J/kg), its derivatives with p and t (output)
    ! irng: flag of out-of-range status (optional)
    !     0: p and t inside the range, 1-9: out-of-range as below
    !
    !     1  |     2     | 3
    !    ---------------------
    !        |           |
    !     4  | t^  0     | 5       9: invalid number (negative or NaN)
    !        |  |        | 
    !        |  +--->p   |
    !    ---------------------
    !     6  |     7     | 8
    !
    */

    if((p <= 0.0 || t <= 0.0) || 
        ((p <= 0.0) && (p > 0.0)) || 
        ((t <= 0.0) && (t > 0.0)))
    {
        printf("wrsteamtab: invalid args range !\n");
        printf("p = %f\n", p);
        printf("t = %f\n", t);

        if (irng != -1) {irng = 9;} // if(present(irng)): irng = 9
        return; // return r, drdp, drdt, h, dhdp, dhdt, irng;
    }
    p1 = log10(p/pcwater);
    t1 = t/tcwater-1.0;

    // ! water
    if(iphase == 0) 
    {
        /* """ // ! density """ */
        iptabl(p1,t1,mp1, mt1l, jmaxl, tabdenl, p1mesh, t1meshl, 
                        f,fx,fy, kp,kt, ifp,ift);
        if (ifp != 0 or ift != 0)
        {
            eptabl(p1,t1,kp,kt,ifp,ift, mp1,mt1l,jmaxl,tabdenl,p1mesh,t1meshl,0, 
                 f,fx,fy);
        }
        // end if       
        r = f;
        drdp = fx/p/log(10);
        drdt = fy/tcwater;

        /* """// ! enthalpy""" */
        iptabl(p1,t1,mp1, mt1l, jmaxl, tabentl, p1mesh, t1meshl,f,fx,fy,kp,kt,ifp,ift);
        if (ifp != 0 || ift != 0)
        {
            eptabl(p1,t1,kp,kt,ifp,ift, mp1,mt1l,jmaxl,tabentl,p1mesh,t1meshl,0, 
                   f,fx,fy);
        }
        // end if        
        h = f;
        dhdp = fx/p/log(10);
        dhdt = fy/tcwater;
    }
   
    // ! steam
    else 
    {
        /* """ // ! density """ */
        iptabg(p1,t1,mp1, mt1g, jmaxg, tabdeng, p1mesh, t1meshg, 
              f,fx,fy, kp,kt, ifp,ift);
        if (ifp != 0 || ift != 0)
        {
            eptabg(p1,t1,kp,kt,ifp,ift, 
                  mp1,mt1g,jmaxg,tabdeng,p1mesh,t1meshg,1,
                  f,fx,fy);
        }
        // end if       
        r = f;
        drdp = fx/p/log(10);
        drdt = fy/tcwater;

        /* """// ! enthalpy""" */
        iptabg(p1,t1,mp1, mt1g, jmaxg, tabentg, 
            p1mesh, t1meshg, f,fx,fy, kp,kt,ifp,ift);
        if (ifp != 0 || ift != 0)
        {
            eptabg(p1,t1,kp,kt,ifp,ift,
                mp1,mt1g,jmaxg,tabentg,p1mesh,t1meshg,2, 
                f,fx,fy);
        }
        // end if        
        h = f;
        dhdp = fx/p/log(10);
        dhdt = fy/tcwater;
    }
    // end if


    // ! check rho<0
    // !  rho=2.17e-3 at 1e3Pa,1000K -> rho=2e-4 at 1e3Pa, 10000K
    // !  let's take 1e-4 as the minimum reasonable value
    if(r < 1.e-4) {
        r = 1.e-4;
        drdp = 0.0;
        drdt = 0.0;
    }
    // end if
    // !
    // ! set out-of-range flag
    // !
    // if(irng == 1) // 1 is present
    // {
        if(ift == 1)
        {
            if(ifp == -1) {irng=1;}
            else if(ifp == 0) {irng=2;}
            else {irng=3;}
            // end if
        }
        else if(ift == -1)
        {
            if(ifp == -1) {irng=6;}
            else if(ifp == 0) {irng=7;}
            else {irng=8;}
            // end if
        }
        else
        {
            if(ifp == -1) {irng=4;}
            else if(ifp == 0) {irng=0;}
            else {irng=5;}
        }
        // end if
    // }
    // end if
    // end if

    return;
}
/*###############################################################################*/ 

/* 
! programs to give tsat, viscosity, thermal conductivity,
! and surface tension
! =======================================================
! function saturation temperature of water
!
! the function revised so that more accurate near 373K,1.013bar
! Tsat=a*p1**3+b*p1**2+c*p1+d+e/p1
! a=4.38027  +/- 0.3118
! b=-57.3072 +/- 6.686
! c=336.053  +/- 52.19
! d=-560.866 +/- 175.3
! e=690.381  +/- 213
*/

/*###############################################################################*/ 
double tsatwaterf(double p)
{
    /* 
    ! deviation from JSME steam table:
    !  Tsat(P) +-0.3K@278-373K +0.5K@453K +1.5K@530K 0@595K -1.5K@637K
    !  Psat(T) +2%@370K +1%@400K -2.5%@500K -2%@550K
    */
    double p1, tsatwater;
    p1 = max(p, 1.0);
    p1 = log10(p1);
    tsatwater = 690.381/p1 - 560.866 + 336.053*p1
                -57.3072*pow(p1,2) + 4.38027*pow(p1,3);
    return tsatwater;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
double psatwaterf(double t)
{
    /* ! psat from t via bi-section method */
    double f0, f1, fa, p0, p1, pa, psatwater;
    double tsatwater; 
    int i; 
    p0 = 614;  // ! min of the available range
    p1 = 25e+6;  // ! initial 1upper lim. around the crit. press
    tsatwater = tsatwaterf(p0);
    f0 = tsatwater - t;
    
    i=1;
    while (i<21) // ! in case p1 is still too small
    {
        tsatwater=tsatwaterf(p1);
        f1 = tsatwater - t;
        if(f1 > 0.0) {break;}
        p1 = p1*2.0;
        i+=1;
    }
    // end do
    
    while (1)
    {
        pa = 0.5*(p0+p1);
        tsatwater=tsatwaterf(pa);
        fa = tsatwater - t;
        if(f0*fa <= 0.0)
        {
            f1 = fa;
            p1 = pa;
        }
        else if(f1*fa <= 0.0)
        {
            f0 = fa;
            p0 = pa;
        }
        else
        {
            printf ("*psatwater* bi-section failed t,ta(K)=%f, %f\n",t, fa+t);
            exit(EXIT_FAILURE);
        }
        // end if
        if(fabs(p1-p0)/max(fabs(p1+p0),1.0) < 1e-6) {break;}
    }
    // end do
    psatwater = pa;
    return psatwater;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
double dtsdpwaterf(double p)
{
    /*
    ! function derivative of saturation temperature of water 
    ! vs pressure
    !  k.moriyama apr,2000
    ! algebraic derivative from tsatwater

    ! p: pressure (Pa)
    ! dtsdpwater: dTsat/dp (K/Pa)    
    */
    double p1; 
    double dtsdpwater;

    p1 = max( p, 1.0);
    p1 = log10(p1);
    if(p <= 1.e8 && p1 > 0.0)
    {
        dtsdpwater =  
            - 690.381/pow(p1,2) + 336.053 -2.0*57.3072*p1 + 3.0*4.38027*pow(p1,2);
        dtsdpwater = dtsdpwater/max(p,1.0)/log(10.0);
    }
    else{
        dtsdpwater = 0.0;
    }
    // endif

    return dtsdpwater;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
double viscwaterf(double roi, double ti)
{
    /* 
    ! viscosity (Internationally Recommended Interpolation Equation
    !            for Dynamic Viscosity (1975))
    !  valid for 273.15<t<1073 K, 0<ro<1050 
    !
    ! k.moriyama 99/9/28
    !
    ! ti:  temperature K
    ! roi: density kg/m^3
    ! viscwater: viscosity Pa s
    */
    int i, j;
    double ro, t, ra, ta, mu, mu0, aaa, bbb;
    double a[4]= {0.0181583, 0.0177624, 0.0105287, -0.0036744};
    double b[6][5] = {
        {0.501938,  0.235622, -0.274637,  0.145831, -0.0270448},
        {0.162888,  0.789393, -0.743539,  0.263129, -0.0253093},
        {-0.130356,  0.673665, -0.959456,  0.347247, -0.0267758},
        {0.907919,  1.207552, -0.687343,  0.213486, -0.0822904},
        {-0.551119,  0.0670665,-0.497089,  0.100754,  0.0602253},
        {0.146543, -0.0843370, 0.195286, -0.032932, -0.0202595}
    };
    double viscwater;
    ra=317.763; 
    ta=647.27;

    ro = roi;
    t = ti;
    // !
    if(t < 273.0){t=273.00;}
    // end if

    if(ro > 1050.0){ro=1050.0;}
    // end if

    aaa=0.0;
    i=0;
    while (i<4){
        aaa=aaa + a[i]*pow((ta/t),i);
        i+=1;
    }
    // end do
    mu0=pow((t/ta),0.5) / aaa;
    bbb=0.0;
    i=0;
    while (i<6){
        j=0;
        while (j<5){
            bbb=bbb + b[i][j]*pow((ta/t-1.0),i)*pow((ro/ra-1.0),j);
            j+=1;
        }
        // end do
        i+=1;
    }
    // end do

    mu = mu0*exp(ro/ra*bbb);
    viscwater=mu*1.e-6;

    return viscwater;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
double tconwaterf (double roi,double ti)
{
    /* 
    ! thermal conductivity 
    ! (Internationally Recommended Interpolation Equation
    !    for Thermal Conductivity for Industrial Use (1977))
    !
    ! ti:  temperature K
    ! roi: density kg/m^3
    ! tconwater: thermal conductivity W/mK
    !
    ! k.moriyama 99/9/28 
    */
    int i;
    double ro, t, ra, ta, lam, lam0, lamb, dlam, qq, rr, ss, dta, aaa;
    double bb[2] = {-1.71587e-1, 2.39219};
    double cc[6] = {6.42857e-1, -4.11717, -6.17937, 3.08976e-3, 8.22994e-2, 1.00932e1};
    double a[4] = {1.02811e-2, 2.99621e-2, 1.56146e-2, -4.22464e-3};
    double b[3] = {-3.97070e-1, 4.00302e-1, 1.06000};
    double d[4] = {7.01309e-2, 1.18520e-2, 1.69937e-3, -1.02000};
    double tconwater;

    ra = 317.7; 
    ta = 647.30;

    ro = roi;
    t = ti;

    // ! give a limiter
    if(t > 2500.0){t = 2500.0;}
    // end if

    // !m030107  t<274 makes invalid value. set a limiter.
    if(t < 274.0){t = 274.0;}
    // end if

    aaa=0.0;
    i=0;
    while (i<4){
        aaa=aaa + a[i]*pow((t/ta),i);
        i+=1;
    }
    // end do

    lam0=pow((t/ta),0.5) * aaa;
    lamb=b[0]+b[1]*(ro/ra)+b[2]*exp(bb[0]*pow((ro/ra+bb[1]),2));
    dta=fabs(t/ta-1.0)+cc[3];
    qq=2.0+cc[4]/pow(dta,0.6);
    rr=qq+1.0;
    if(t/ta >= 1.0){
        ss=1.0/dta;
    }
    else{
        ss=cc[5]/pow(dta,0.6);
    }   
    // end if
    dlam=(d[0]*pow((ta/t),10)+d[1])*pow((ro/ra),1.8)
        *exp(cc[0]*(1.0-pow((ro/ra),2.8))) 
        +d[2]*ss*pow((ro/ra),qq)*exp(qq/rr*(1.0-pow((ro/ra),rr))) 
        +d[3]*exp(cc[1]*pow((t/ta),1.5)+cc[2]*pow((ra/ro),5));
    lam=lam0+lamb+dlam;
    tconwater=lam;

    return tconwater;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
double stenwaterf(double t)
{
    /* 
    ! surface tension (Internationally Recommended Interpolation equation
    !    for Surface Tension (1975))
    ! t:  temperature K
    ! stenwater: surface tension N/m
    ! k.moriyama 99/9/28
    ! m30901 give min. 5mN/m to avoid numeric errors.
    !        it is about the value at 640K. 
    */    
    double const tc = 647.15; 
    double stenwater;

    if(t >= tc){
    // !m030901         stenwater=0.d0
        stenwater=5.e-3;
        return stenwater;
    }
    // end if
    stenwater = 0.2358*pow((1.0-t/tc),1.256)*(1.0-0.625*(1.0-t/tc));
    // ! m030901
    stenwater = max(stenwater, 5.e-3);
    
    return stenwater;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
double hfgwaterf(double p)
{
    /*     
    ! latent heat of evaporation
    !  note: this is made from the result of wrsteam calculation
    !        for usage in heat ransfer correlation etc.
    !        consistency with the wrsteamtab is not perfect.
    !        but somehow this is more precise with real physics; 
    !        i.e. upper limited for p<1e3 and drops to 0 at pcrit. 
    */
    double const pc = 22.12e6, // crit. pressure
    x1 = 0.467576500991882, a = 7576210.54878548, b = -22923268.2138487, 
    c = 35469825.3810538, d = -18542657.5515397, pa = -0.05, ha = 835952.126288328,
    pb = -4.5, hb = 2508645.26205592;
    double p1; 
    double hfgwater;

    p1 = log10(p/pc);
    if(p1 < pb){
        hfgwater = hb;
    }
    else if(p1 < pa){
        hfgwater = a*pow((x1-p1),0.7) + b*pow((x1-p1),0.5) + c*pow((x1-p1),0.2) + d;
    }
    else if(p1 < 0.0){
        hfgwater = ha/pa*p1;
    }
    else{
        hfgwater = 0.0;
    }
    // end if

    return hfgwater;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
double betawaterf(int n, double t)
{
    /* 
    ! coefficient of cubic expansion
    !  note: this is made separately from calc. by steamsi package
    !        for usage in heat ransfer correlation etc.
    !        (espacially for natural convection heat transfer)
    !        consistency with the wrsteamtab is totally not maintained.
    !        but this is more precise with real physics. 
    */
    double const a=6.140e-6, b=-1.559e-3, t1min=280.0, t1max=550.0;
    double t1;
    double betawater;

    if(n == 0){ // ! water
        t1 = min( max(t, t1min), t1max );
        betawater = a*t1+b;
    }
    else{              // ! steam (ideal gas assumption)
        t1 = max(t, t1min);
        betawater = 1.0/t1;
    }
    // end if

    return betawater;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void wrsttabset()
{
    /* 
    ! definitions requiered for table handling ------------------
    !
    ! -----------------------------------------------------
    ! table data was made by mkeos_water.f program by k.moriyama
    ! referring steamsi.f program.
    ! steamsi.f is a modified version of steam.f program
    ! originally written by k.kobayashi (1976) referring
    ! JSME steam tables (1968) (C) Japan Association of Mechanical 
    ! Engineers.
    !
    ! out.prog generated by mkeos_water.f comes here 
    */

    //  ! pressure mesh (gas/liq)
    p1mesh[ 0] = -4.344790E+000;
    p1mesh[ 1] = -3.500000E+000;
    p1mesh[ 2] = -3.000000E+000;
    p1mesh[ 3] = -2.500000E+000;
    p1mesh[ 4] = -2.000000E+000;
    p1mesh[ 5] = -1.500000E+000;
    p1mesh[ 6] = -1.000000E+000;
    p1mesh[ 7] = -6.000000E-001;
    p1mesh[ 8] = -2.000000E-001;
    p1mesh[ 9] =  2.000000E-001;
    p1mesh[10] =  6.552140E-001;

    //  !
    //  ! gas temperature mesh limits
    //  !
    jmaxg[ 0] = 12;
    jmaxg[ 1] = 12;
    jmaxg[ 2] = 12;
    jmaxg[ 3] = 11;
    jmaxg[ 4] = 10;
    jmaxg[ 5] =  9;
    jmaxg[ 6] =  8;
    jmaxg[ 7] =  7;
    jmaxg[ 8] =  6;
    jmaxg[ 9] =  5;
    jmaxg[10] =  4;
    //  !
    //  ! liquid temperature mesh limits
    //  !
    jmaxl[ 0] =  2;
    jmaxl[ 1] =  3;
    jmaxl[ 2] =  4;
    jmaxl[ 3] =  5;
    jmaxl[ 4] =  6;
    jmaxl[ 5] =  7;
    jmaxl[ 6] =  8;
    jmaxl[ 7] =  9;
    jmaxl[ 8] = 10;
    jmaxl[ 9] = 10;
    jmaxl[10] = 10;
    //  !
    //  ! gas temperature mesh
    //  !
    t1meshg[ 0] =  5.448000E-001;
    t1meshg[ 1] =  4.000000E-001;
    t1meshg[ 2] =  2.173000E-001;
    t1meshg[ 3] =  1.000000E-001;
    t1meshg[ 4] =  0.000000E+000;
    t1meshg[ 5] = -1.400000E-001;
    t1meshg[ 6] = -2.398000E-001;
    t1meshg[ 7] = -3.000000E-001;
    t1meshg[ 8] = -3.866000E-001;
    t1meshg[ 9] = -4.200000E-001;
    t1meshg[10] = -4.790000E-001;
    t1meshg[11] = -5.000000E-001;
    t1meshg[12] = -5.658000E-001;
    //  !
    //  ! liquid temperature mesh
    //  !
    t1meshl[ 0] = -5.767000E-001;
    t1meshl[ 1] = -5.300000E-001;
    t1meshl[ 2] = -4.850000E-001;
    t1meshl[ 3] = -4.600000E-001;
    t1meshl[ 4] = -3.898000E-001;
    t1meshl[ 5] = -3.400000E-001;
    t1meshl[ 6] = -2.431000E-001;
    t1meshl[ 7] = -1.800000E-001;
    t1meshl[ 8] = -1.000000E-001;
    t1meshl[ 9] =  3.000000E-002;
    t1meshl[10] =  2.173000E-001;
    

;
    // !;
    // ! gas table [density];
    // !;

    tabdeng[ 0][ 0][ 0] =  2.16689241609823455073E-003;
    tabdeng[ 0][ 1][ 0] =  2.39101451560083970091E-003;
    tabdeng[ 0][ 2][ 0] =  2.74988175156841430091E-003;
    tabdeng[ 0][ 3][ 0] =  3.04312969121501696729E-003;
    tabdeng[ 0][ 4][ 0] =  3.34745789687756686859E-003;
    tabdeng[ 0][ 5][ 0] =  3.89243755623770509508E-003;
    tabdeng[ 0][ 6][ 0] =  4.40351520672352045316E-003;
    tabdeng[ 0][ 7][ 0] =  4.78230394867777599405E-003;
    tabdeng[ 0][ 8][ 0] =  5.45774425325034392292E-003;
    tabdeng[ 0][ 9][ 0] =  5.77222888977193485799E-003;
    tabdeng[ 0][10][ 0] =  6.42652155630545757203E-003;
    tabdeng[ 0][11][ 0] =  6.69678781894629104182E-003;
    tabdeng[ 0][12][ 0] =  7.71377542790852716026E-003;
    tabdeng[ 1][ 0][ 0] =  1.51577406138793738255E-002;
    tabdeng[ 1][ 1][ 0] =  1.67256570163615107816E-002;
    tabdeng[ 1][ 2][ 0] =  1.92363797744052981986E-002;
    tabdeng[ 1][ 3][ 0] =  2.12881936128708620926E-002;
    tabdeng[ 1][ 4][ 0] =  2.34177586854883885981E-002;
    tabdeng[ 1][ 5][ 0] =  2.72321380484774229480E-002;
    tabdeng[ 1][ 6][ 0] =  3.08108133466362087149E-002;
    tabdeng[ 1][ 7][ 0] =  3.34647812075495112993E-002;
    tabdeng[ 1][ 8][ 0] =  3.82027608031849444381E-002;
    tabdeng[ 1][ 9][ 0] =  4.04122430650496744509E-002;
    tabdeng[ 1][10][ 0] =  4.50193962354188195740E-002;
    tabdeng[ 1][11][ 0] =  4.69275938761514466169E-002;
    tabdeng[ 1][12][ 0] =  5.16425576502705444004E-002;
    tabdeng[ 2][ 0][ 0] =  4.79349725721072983387E-002;
    tabdeng[ 2][ 1][ 0] =  5.28945491515248542025E-002;
    tabdeng[ 2][ 2][ 0] =  6.08376099321372154627E-002;
    tabdeng[ 2][ 3][ 0] =  6.73302757048018646335E-002;
    tabdeng[ 2][ 4][ 0] =  7.40707649266291318080E-002;
    tabdeng[ 2][ 5][ 0] =  8.61506877190813535883E-002;
    tabdeng[ 2][ 6][ 0] =  9.74968327296403536320E-002;
    tabdeng[ 2][ 7][ 0] =  1.05924159286228378174E-001;
    tabdeng[ 2][ 8][ 0] =  1.21013834256805868916E-001;
    tabdeng[ 2][ 9][ 0] =  1.28079016560593117458E-001;
    tabdeng[ 2][10][ 0] =  1.42895432292380053507E-001;
    tabdeng[ 2][11][ 0] =  1.46419123184130889337E-001;
    tabdeng[ 2][12][ 0] =  1.47219741879111132299E-001;
    tabdeng[ 3][ 0][ 0] =  1.51603573773206795616E-001;
    tabdeng[ 3][ 1][ 0] =  1.67301032729939236354E-001;
    tabdeng[ 3][ 2][ 0] =  1.92453803446372279096E-001;
    tabdeng[ 3][ 3][ 0] =  2.13028076634167257142E-001;
    tabdeng[ 3][ 4][ 0] =  2.34405645843175608478E-001;
    tabdeng[ 3][ 5][ 0] =  2.72784594239217470513E-001;
    tabdeng[ 3][ 6][ 0] =  3.08961457842393039108E-001;
    tabdeng[ 3][ 7][ 0] =  3.35965247895937180189E-001;
    tabdeng[ 3][ 8][ 0] =  3.84787200357978420584E-001;
    tabdeng[ 3][ 9][ 0] =  4.07947939963371586369E-001;
    tabdeng[ 3][10][ 0] =  4.44384068151180311457E-001;
    tabdeng[ 3][11][ 0] =  4.46626756665510826760E-001;
    tabdeng[ 4][ 0][ 0] =  4.79611414980400774244E-001;
    tabdeng[ 4][ 1][ 0] =  5.29390307307490703970E-001;
    tabdeng[ 4][ 2][ 0] =  6.09277195311073160866E-001;
    tabdeng[ 4][ 3][ 0] =  6.74767123754525965929E-001;
    tabdeng[ 4][ 4][ 0] =  7.42996189315748178927E-001;
    tabdeng[ 4][ 5][ 0] =  8.66180347713465526027E-001;
    tabdeng[ 4][ 6][ 0] =  9.83668215238930554456E-001;
    tabdeng[ 4][ 7][ 0] =  1.07283565707626249441E+000;
    tabdeng[ 4][ 8][ 0] =  1.23949584536434387338E+000;
    tabdeng[ 4][ 9][ 0] =  1.29419393390306591840E+000;
    tabdeng[ 4][10][ 0] =  1.32902393256168172009E+000;
    tabdeng[ 5][ 0][ 0] =  1.51865305913400772120E+000;
    tabdeng[ 5][ 1][ 0] =  1.67746461503851751473E+000;
    tabdeng[ 5][ 2][ 0] =  1.93358291474575061386E+000;
    tabdeng[ 5][ 3][ 0] =  2.14502304779118269096E+000;
    tabdeng[ 5][ 4][ 0] =  2.36721346012622957389E+000;
    tabdeng[ 5][ 5][ 0] =  2.77605307262863476581E+000;
    tabdeng[ 5][ 6][ 0] =  3.18271329116733170395E+000;
    tabdeng[ 5][ 7][ 0] =  3.51125223264085573049E+000;
    tabdeng[ 5][ 8][ 0] =  3.98340713767065901507E+000;
    tabdeng[ 5][ 9][ 0] =  4.06221664786154779136E+000;
    tabdeng[ 6][ 0][ 0] =  4.82230290376491588233E+000;
    tabdeng[ 6][ 1][ 0] =  5.33865132495352057163E+000;
    tabdeng[ 6][ 2][ 0] =  6.18439813896842416341E+000;
    tabdeng[ 6][ 3][ 0] =  6.89869227476474744520E+000;
    tabdeng[ 6][ 4][ 0] =  7.67201124339578122147E+000;
    tabdeng[ 6][ 5][ 0] =  9.20487811904215291747E+000;
    tabdeng[ 6][ 6][ 0] =  1.10363885811016935179E+001;
    tabdeng[ 6][ 7][ 0] =  1.22884199507020674957E+001;
    tabdeng[ 6][ 8][ 0] =  1.34605635508878691553E+001;
    tabdeng[ 7][ 0][ 0] =  1.22238569667049485901E+001;
    tabdeng[ 7][ 1][ 0] =  1.36020600655522621025E+001;
    tabdeng[ 7][ 2][ 0] =  1.59404756869264545571E+001;
    tabdeng[ 7][ 3][ 0] =  1.80298004009703696227E+001;
    tabdeng[ 7][ 4][ 0] =  2.04908325826366066735E+001;
    tabdeng[ 7][ 5][ 0] =  2.67508731427973849293E+001;
    tabdeng[ 7][ 6][ 0] =  3.34358404188304234594E+001;
    tabdeng[ 7][ 7][ 0] =  3.64214799325265730090E+001;
    tabdeng[ 8][ 0][ 0] =  3.14094629591204146379E+001;
    tabdeng[ 8][ 1][ 0] =  3.54338048309352871001E+001;
    tabdeng[ 8][ 2][ 0] =  4.29850163555816138228E+001;
    tabdeng[ 8][ 3][ 0] =  5.11169156740205750111E+001;
    tabdeng[ 8][ 4][ 0] =  6.45730519149321509076E+001;
    tabdeng[ 8][ 5][ 0] =  9.23241928890473531055E+001;
    tabdeng[ 8][ 6][ 0] =  1.08636391116324531936E+002;
    tabdeng[ 9][ 0][ 0] =  8.34245957425389121909E+001;
    tabdeng[ 9][ 1][ 0] =  9.80315242377574804777E+001;
    tabdeng[ 9][ 2][ 0] =  1.35808467175770118729E+002;
    tabdeng[ 9][ 3][ 0] =  2.11582561751117424365E+002;
    tabdeng[ 9][ 4][ 0] =  3.70162070630437824548E+002;
    tabdeng[ 9][ 5][ 0] =  4.51756184820851672157E+002;
    tabdeng[10][ 0][ 0] =  2.65808659214430804241E+002;
    tabdeng[10][ 1][ 0] =  3.37763153098218595005E+002;
    tabdeng[10][ 2][ 0] =  4.82861515601031101141E+002;
    tabdeng[10][ 3][ 0] =  6.53821683291788758652E+002;
    tabdeng[10][ 4][ 0] =  8.61112163947147564613E+002;
    tabdeng[ 0][ 0][ 1] =  4.98946778854920426627E-003;
    tabdeng[ 0][ 1][ 1] =  5.50553754846578105586E-003;
    tabdeng[ 0][ 2][ 1] =  6.33188367833453259748E-003;
    tabdeng[ 0][ 3][ 1] =  7.00714131030832909630E-003;
    tabdeng[ 0][ 4][ 1] =  7.70792560960817559057E-003;
    tabdeng[ 0][ 5][ 1] =  8.96290984077656276152E-003;
    tabdeng[ 0][ 6][ 1] =  1.01399108166476056847E-002;
    tabdeng[ 0][ 7][ 1] =  1.10123412586489288073E-002;
    tabdeng[ 0][ 8][ 1] =  1.25683253095325853099E-002;
    tabdeng[ 0][ 9][ 1] =  1.32929805416580962274E-002;
    tabdeng[ 0][10][ 1] =  1.48012016992796653636E-002;
    tabdeng[ 0][11][ 1] =  1.54244759442248944442E-002;
    tabdeng[ 0][12][ 1] =  1.77717935994064490934E-002;
    tabdeng[ 1][ 0][ 1] =  3.49026565672187968903E-002;
    tabdeng[ 1][ 1][ 1] =  3.85133853448873866854E-002;
    tabdeng[ 1][ 2][ 1] =  4.42957024366094076484E-002;
    tabdeng[ 1][ 3][ 1] =  4.90216126764552381778E-002;
    tabdeng[ 1][ 4][ 1] =  5.39272085536526210592E-002;
    tabdeng[ 1][ 5][ 1] =  6.27161269932050291498E-002;
    tabdeng[ 1][ 6][ 1] =  7.09661979471629678073E-002;
    tabdeng[ 1][ 7][ 1] =  7.70888305546548607827E-002;
    tabdeng[ 1][ 8][ 1] =  8.80341397071871117896E-002;
    tabdeng[ 1][ 9][ 1] =  9.31477016793497425429E-002;
    tabdeng[ 1][10][ 1] =  1.03838097728375794437E-001;
    tabdeng[ 1][11][ 1] =  1.08279732986442167597E-001;
    tabdeng[ 1][12][ 1] =  8.62275036616158285785E-002;
    tabdeng[ 2][ 0][ 1] =  1.10381046491237838625E-001;
    tabdeng[ 2][ 1][ 1] =  1.21805573172319964170E-001;
    tabdeng[ 2][ 2][ 1] =  1.40106793921171646211E-001;
    tabdeng[ 2][ 3][ 1] =  1.55071061321960385060E-001;
    tabdeng[ 2][ 4][ 1] =  1.70612544948872985051E-001;
    tabdeng[ 2][ 5][ 1] =  1.98487596087270679046E-001;
    tabdeng[ 2][ 6][ 1] =  2.24712267837526818992E-001;
    tabdeng[ 2][ 7][ 1] =  2.44234442685705432918E-001;
    tabdeng[ 2][ 8][ 1] =  2.79342391172537951594E-001;
    tabdeng[ 2][ 9][ 1] =  2.95876765285877252332E-001;
    tabdeng[ 2][10][ 1] =  3.30838292459378346955E-001;
    tabdeng[ 2][11][ 1] =  2.89686384245475603283E-001;
    tabdeng[ 2][12][ 1] =  2.96081233253746523015E-001;
    tabdeng[ 3][ 0][ 1] =  3.49147073439048438726E-001;
    tabdeng[ 3][ 1][ 1] =  3.85338631282030774639E-001;
    tabdeng[ 3][ 2][ 1] =  4.43371625709419403538E-001;
    tabdeng[ 3][ 3][ 1] =  4.90889439985371323072E-001;
    tabdeng[ 3][ 4][ 1] =  5.40323165908118463463E-001;
    tabdeng[ 3][ 5][ 1] =  6.29298709278521939403E-001;
    tabdeng[ 3][ 6][ 1] =  7.13608710294686976461E-001;
    tabdeng[ 3][ 7][ 1] =  7.76997956733024852660E-001;
    tabdeng[ 3][ 8][ 1] =  8.93226498604199559139E-001;
    tabdeng[ 3][ 9][ 1] =  9.49410292179885928476E-001;
    tabdeng[ 3][10][ 1] =  8.75116250975822684843E-001;
    tabdeng[ 3][11][ 1] =  9.11144149680044090900E-001;
    tabdeng[ 4][ 0][ 1] =  1.10501563248749934409E+000;
    tabdeng[ 4][ 1][ 1] =  1.22010482738943792924E+000;
    tabdeng[ 4][ 2][ 1] =  1.40522116559315324302E+000;
    tabdeng[ 4][ 3][ 1] =  1.55746437212771993153E+000;
    tabdeng[ 4][ 4][ 1] =  1.71669186235978954436E+000;
    tabdeng[ 4][ 5][ 1] =  2.00654166023288516030E+000;
    tabdeng[ 4][ 6][ 1] =  2.28777064217562564608E+000;
    tabdeng[ 4][ 7][ 1] =  2.50642463205411392124E+000;
    tabdeng[ 4][ 8][ 1] =  2.93489395383507734039E+000;
    tabdeng[ 4][ 9][ 1] =  2.59557368357889117760E+000;
    tabdeng[ 4][10][ 1] =  2.66344320666618283866E+000;
    tabdeng[ 5][ 0][ 1] =  3.50352544044936564660E+000;
    tabdeng[ 5][ 1][ 1] =  3.87392000648424561504E+000;
    tabdeng[ 5][ 2][ 1] =  4.47548640767860828049E+000;
    tabdeng[ 5][ 3][ 1] =  4.97712975739440111056E+000;
    tabdeng[ 5][ 4][ 1] =  5.51083883778590610802E+000;
    tabdeng[ 5][ 5][ 1] =  6.52035125375197566200E+000;
    tabdeng[ 5][ 6][ 1] =  7.58736220089817958723E+000;
    tabdeng[ 5][ 7][ 1] =  8.52670584823549937425E+000;
    tabdeng[ 5][ 8][ 1] =  8.04075121539018411454E+000;
    tabdeng[ 5][ 9][ 1] =  8.47651717225503631425E+000;
    tabdeng[ 6][ 0][ 1] =  1.11708172429429328787E+001;
    tabdeng[ 6][ 1][ 1] =  1.24078421737546449322E+001;
    tabdeng[ 6][ 2][ 1] =  1.44783942368540738954E+001;
    tabdeng[ 6][ 3][ 1] =  1.62834289277879022961E+001;
    tabdeng[ 6][ 4][ 1] =  1.83214454815416125655E+001;
    tabdeng[ 6][ 5][ 1] =  2.28054829685534521388E+001;
    tabdeng[ 6][ 6][ 1] =  2.95215290117507400680E+001;
    tabdeng[ 6][ 7][ 1] =  2.65819650240093494631E+001;
    tabdeng[ 6][ 8][ 1] =  2.98678744374786582227E+001;
    tabdeng[ 7][ 0][ 1] =  2.85710630365191384783E+001;
    tabdeng[ 7][ 1][ 1] =  3.20627667944593284233E+001;
    tabdeng[ 7][ 2][ 1] =  3.83120034835747986790E+001;
    tabdeng[ 7][ 3][ 1] =  4.43860645773335775743E+001;
    tabdeng[ 7][ 4][ 1] =  5.24557146250910122376E+001;
    tabdeng[ 7][ 5][ 1] =  8.07211244453653051778E+001;
    tabdeng[ 7][ 6][ 1] =  8.24757301768929096397E+001;
    tabdeng[ 7][ 7][ 1] =  9.40833348851131745505E+001;
    tabdeng[ 8][ 0][ 1] =  7.50359515184042749070E+001;
    tabdeng[ 8][ 1][ 1] =  8.66104714631855046036E+001;
    tabdeng[ 8][ 2][ 1] =  1.11522661486184347268E+002;
    tabdeng[ 8][ 3][ 1] =  1.45122563971142483297E+002;
    tabdeng[ 8][ 4][ 1] =  2.26302920458248507884E+002;
    tabdeng[ 8][ 5][ 1] =  2.47145474285884574783E+002;
    tabdeng[ 8][ 6][ 1] =  2.93527023310577646953E+002;
    tabdeng[ 9][ 0][ 1] =  2.09423111579136900673E+002;
    tabdeng[ 9][ 1][ 1] =  2.62569849864185584920E+002;
    tabdeng[ 9][ 2][ 1] =  4.50095480535389469878E+002;
    tabdeng[ 9][ 3][ 1] =  1.34410874604856371661E+003;
    tabdeng[ 9][ 4][ 1] =  1.30164217311927973242E+003;
    tabdeng[ 9][ 5][ 1] =  1.55001448537313694942E+003;
    tabdeng[10][ 0][ 1] =  6.49498487912990526638E+002;
    tabdeng[10][ 1][ 1] =  7.75512645359187899885E+002;
    tabdeng[10][ 2][ 1] =  7.01958186153046312938E+002;
    tabdeng[10][ 3][ 1] =  5.98885582790692183153E+002;
    tabdeng[10][ 4][ 1] =  8.55365654985638684593E+002;
    tabdeng[ 0][ 0][ 2] = -1.40271844528354143987E-003;
    tabdeng[ 0][ 1][ 2] = -1.70789801256991430507E-003;
    tabdeng[ 0][ 2][ 2] = -2.25906583809502197260E-003;
    tabdeng[ 0][ 3][ 2] = -2.76659314750605372923E-003;
    tabdeng[ 0][ 4][ 2] = -3.34764610663439588639E-003;
    tabdeng[ 0][ 5][ 2] = -4.52655010703661986327E-003;
    tabdeng[ 0][ 6][ 2] = -5.79360559331902613689E-003;
    tabdeng[ 0][ 7][ 2] = -6.83368520609645754144E-003;
    tabdeng[ 0][ 8][ 2] = -8.90215072932821344098E-003;
    tabdeng[ 0][ 9][ 2] = -9.95896237060274990005E-003;
    tabdeng[ 0][10][ 2] = -1.23492711121849750577E-002;
    tabdeng[ 0][11][ 2] = -1.34124037375139299605E-002;
    tabdeng[ 0][12][ 2] = -1.78121529133067613149E-002;
    tabdeng[ 1][ 0][ 2] = -9.81296133843525923701E-003;
    tabdeng[ 1][ 1][ 2] = -1.19483903955406173258E-002;
    tabdeng[ 1][ 2][ 2] = -1.58056766910615320576E-002;
    tabdeng[ 1][ 3][ 2] = -1.93583653817584674939E-002;
    tabdeng[ 1][ 4][ 2] = -2.34269707709030446108E-002;
    tabdeng[ 1][ 5][ 2] = -3.16877994205399873828E-002;
    tabdeng[ 1][ 6][ 2] = -4.05804076430420837540E-002;
    tabdeng[ 1][ 7][ 2] = -4.78962200268181759188E-002;
    tabdeng[ 1][ 8][ 2] = -6.25074484348640962983E-002;
    tabdeng[ 1][ 9][ 2] = -7.00128990819822438763E-002;
    tabdeng[ 1][10][ 2] = -8.71145014073218398876E-002;
    tabdeng[ 1][11][ 2] = -9.47845630606958117204E-002;
    tabdeng[ 1][12][ 2] = -4.85274057575139153298E-002;
    tabdeng[ 2][ 0][ 2] = -3.10384471049184214320E-002;
    tabdeng[ 2][ 1][ 2] = -3.77967479983234780350E-002;
    tabdeng[ 2][ 2][ 2] = -5.00093029231758953723E-002;
    tabdeng[ 2][ 3][ 2] = -6.12639867735142173988E-002;
    tabdeng[ 2][ 4][ 2] = -7.41629587780103100014E-002;
    tabdeng[ 2][ 5][ 2] = -1.00400907108556591840E-001;
    tabdeng[ 2][ 6][ 2] = -1.28758940397935506272E-001;
    tabdeng[ 2][ 7][ 2] = -1.52219662703170666163E-001;
    tabdeng[ 2][ 8][ 2] = -1.99579331914889718247E-001;
    tabdeng[ 2][ 9][ 2] = -2.24236976022758854876E-001;
    tabdeng[ 2][10][ 2] = -2.81458119694295516577E-001;
    tabdeng[ 2][11][ 2] = -5.41314890438790163874E-002;
    tabdeng[ 2][12][ 2] = -1.21674573705204160101E-002;
    tabdeng[ 3][ 0][ 2] = -9.82235935087393541298E-002;
    tabdeng[ 3][ 1][ 2] = -1.19650069911034245829E-001;
    tabdeng[ 3][ 2][ 2] = -1.58417204650016185674E-001;
    tabdeng[ 3][ 3][ 2] = -1.94209113397200372475E-001;
    tabdeng[ 3][ 4][ 2] = -2.35329921145879111810E-001;
    tabdeng[ 3][ 5][ 2] = -3.19463168651505091677E-001;
    tabdeng[ 3][ 6][ 2] = -4.11562645509082947193E-001;
    tabdeng[ 3][ 7][ 2] = -4.89124216105071007998E-001;
    tabdeng[ 3][ 8][ 2] = -6.51059676006663812231E-001;
    tabdeng[ 3][ 9][ 2] = -7.38943733450023954035E-001;
    tabdeng[ 3][10][ 2] = -4.96179255967221011225E-001;
    tabdeng[ 3][11][ 2] = -1.06794691158595869696E-001;
    tabdeng[ 4][ 0][ 2] = -3.11324928107844733915E-001;
    tabdeng[ 4][ 1][ 2] = -3.79631182110770126759E-001;
    tabdeng[ 4][ 2][ 2] = -5.03706478435536286753E-001;
    tabdeng[ 4][ 3][ 2] = -6.18922008939958967488E-001;
    tabdeng[ 4][ 4][ 2] = -7.52315649815511355847E-001;
    tabdeng[ 4][ 5][ 2] = -1.03040077983846689591E+000;
    tabdeng[ 4][ 6][ 2] = -1.34764886104026304636E+000;
    tabdeng[ 4][ 7][ 2] = -1.63046374956878414864E+000;
    tabdeng[ 4][ 8][ 2] = -2.28617165156866031239E+000;
    tabdeng[ 4][ 9][ 2] = -9.89162991468589325450E-001;
    tabdeng[ 4][10][ 2] = -1.91514929162454827782E-001;
    tabdeng[ 5][ 0][ 2] = -9.91661507279333998000E-001;
    tabdeng[ 5][ 1][ 2] = -1.21320369243238879164E+000;
    tabdeng[ 5][ 2][ 2] = -1.62060685975004292736E+000;
    tabdeng[ 5][ 3][ 2] = -2.00585112723742797769E+000;
    tabdeng[ 5][ 4][ 2] = -2.46310584603199300702E+000;
    tabdeng[ 5][ 5][ 2] = -3.47820412120348398943E+000;
    tabdeng[ 5][ 6][ 2] = -4.80845683966752268645E+000;
    tabdeng[ 5][ 7][ 2] = -6.22791267227120037830E+000;
    tabdeng[ 5][ 8][ 2] = -4.67635765174273121403E+000;
    tabdeng[ 5][ 9][ 2] = -4.27746950170782014311E-002;
    tabdeng[ 6][ 0][ 2] = -3.20819491517386179424E+000;
    tabdeng[ 6][ 1][ 2] = -3.96556411061288827113E+000;
    tabdeng[ 6][ 2][ 2] = -5.41222529047395006785E+000;
    tabdeng[ 6][ 3][ 2] = -6.86292397354962702138E+000;
    tabdeng[ 6][ 4][ 2] = -8.74089365031807830064E+000;
    tabdeng[ 6][ 5][ 2] = -1.39814048900660861108E+001;
    tabdeng[ 6][ 6][ 2] = -2.47441955030167690666E+001;
    tabdeng[ 6][ 7][ 2] = -1.68515310617797240411E+001;
    tabdeng[ 6][ 8][ 2] = -1.02187599355828986347E+001;
    tabdeng[ 7][ 0][ 2] = -8.46942476007012068351E+000;
    tabdeng[ 7][ 1][ 2] = -1.07097593662347527754E+001;
    tabdeng[ 7][ 2][ 2] = -1.53666602822353528524E+001;
    tabdeng[ 7][ 3][ 2] = -2.07659695241501864871E+001;
    tabdeng[ 7][ 4][ 2] = -2.94759219094044873088E+001;
    tabdeng[ 7][ 5][ 2] = -7.25689370729122629200E+001;
    tabdeng[ 7][ 6][ 2] = -6.13983430079101566434E+001;
    tabdeng[ 7][ 7][ 2] = -3.77923385102343871722E+001;
    tabdeng[ 8][ 0][ 2] = -2.40335539376488540597E+001;
    tabdeng[ 8][ 1][ 2] = -3.22531206429236760869E+001;
    tabdeng[ 8][ 2][ 2] = -5.38311149680566387588E+001;
    tabdeng[ 8][ 3][ 2] = -9.09859730611666464029E+001;
    tabdeng[ 8][ 4][ 2] = -2.10735544539959704480E+002;
    tabdeng[ 8][ 5][ 2] = -1.85709326518828845565E+002;
    tabdeng[ 8][ 6][ 2] = -1.41188433546846084710E+002;
    tabdeng[ 9][ 0][ 2] = -8.07038777155869695434E+001;
    tabdeng[ 9][ 1][ 2] = -1.27384624389230936004E+002;
    tabdeng[ 9][ 2][ 2] = -3.52676098389535184197E+002;
    tabdeng[ 9][ 3][ 2] = -2.12346471512371908830E+003;
    tabdeng[ 9][ 4][ 2] = -1.04812546246268880168E+003;
    tabdeng[ 9][ 5][ 2] = -1.17504740257508956347E+002;
    tabdeng[10][ 0][ 2] = -3.80818735595397527049E+002;
    tabdeng[10][ 1][ 2] = -6.37290041623056822573E+002;
    tabdeng[10][ 2][ 2] = -1.09595621887325864918E+003;
    tabdeng[10][ 3][ 2] = -1.81896565138688924890E+003;
    tabdeng[10][ 4][ 2] = -2.32684396172028664296E+003;
    
    // !;
    // ! gas table [enthalpy];
    // !;
    tabentg[ 0][ 0][ 0] =  3.98911502493808884174E+006;
    tabentg[ 0][ 1][ 0] =  3.77800201012914022431E+006;
    tabentg[ 0][ 2][ 0] =  3.52020123457061825320E+006;
    tabentg[ 0][ 3][ 0] =  3.35976983619108237326E+006;
    tabentg[ 0][ 4][ 0] =  3.22608034500470338389E+006;
    tabentg[ 0][ 5][ 0] =  3.04347720087749650702E+006;
    tabentg[ 0][ 6][ 0] =  2.91635051437951484695E+006;
    tabentg[ 0][ 7][ 0] =  2.84079562208400666714E+006;
    tabentg[ 0][ 8][ 0] =  2.73346624082268262282E+006;
    tabentg[ 0][ 9][ 0] =  2.69246474261438520625E+006;
    tabentg[ 0][10][ 0] =  2.62051796037613740191E+006;
    tabentg[ 0][11][ 0] =  2.59504618638897268102E+006;
    tabentg[ 0][12][ 0] =  2.51564088396537536755E+006;
    tabentg[ 1][ 0][ 0] =  3.98908053919636690989E+006;
    tabentg[ 1][ 1][ 0] =  3.77795723950652917847E+006;
    tabentg[ 1][ 2][ 0] =  3.52013854102862020954E+006;
    tabentg[ 1][ 3][ 0] =  3.35969037151530245319E+006;
    tabentg[ 1][ 4][ 0] =  3.22597966367555223405E+006;
    tabentg[ 1][ 5][ 0] =  3.04332068532766215503E+006;
    tabentg[ 1][ 6][ 0] =  2.91610811240983987227E+006;
    tabentg[ 1][ 7][ 0] =  2.84046100320000248030E+006;
    tabentg[ 1][ 8][ 0] =  2.73289493794284714386E+006;
    tabentg[ 1][ 9][ 0] =  2.69174943680900335312E+006;
    tabentg[ 1][10][ 0] =  2.61943438166987383738E+006;
    tabentg[ 1][11][ 0] =  2.59378469580838596448E+006;
    tabentg[ 1][12][ 0] =  2.51354899455463979393E+006;
    tabentg[ 2][ 0][ 0] =  3.98899353362052282318E+006;
    tabentg[ 2][ 1][ 0] =  3.77784428520161518827E+006;
    tabentg[ 2][ 2][ 0] =  3.51998036459537921473E+006;
    tabentg[ 2][ 3][ 0] =  3.35948986973674083129E+006;
    tabentg[ 2][ 4][ 0] =  3.22572559177025780082E+006;
    tabentg[ 2][ 5][ 0] =  3.04292546342627191916E+006;
    tabentg[ 2][ 6][ 0] =  2.91549534774676058441E+006;
    tabentg[ 2][ 7][ 0] =  2.83961423985024448484E+006;
    tabentg[ 2][ 8][ 0] =  2.73144617113757878542E+006;
    tabentg[ 2][ 9][ 0] =  2.68993362098716944456E+006;
    tabentg[ 2][10][ 0] =  2.61667778841788554564E+006;
    tabentg[ 2][11][ 0] =  2.59066409156172303483E+006;
    tabentg[ 2][12][ 0] =  2.50949269026326062158E+006;
    tabentg[ 3][ 0][ 0] =  3.98871840176267875358E+006;
    tabentg[ 3][ 1][ 0] =  3.77748709210356278345E+006;
    tabentg[ 3][ 2][ 0] =  3.51948013384968973696E+006;
    tabentg[ 3][ 3][ 0] =  3.35885567422259924933E+006;
    tabentg[ 3][ 4][ 0] =  3.22492159488514671102E+006;
    tabentg[ 3][ 5][ 0] =  3.04167240051909070462E+006;
    tabentg[ 3][ 6][ 0] =  2.91354616086293710396E+006;
    tabentg[ 3][ 7][ 0] =  2.83691225476389285177E+006;
    tabentg[ 3][ 8][ 0] =  2.72679395646882755682E+006;
    tabentg[ 3][ 9][ 0] =  2.68408494477462256327E+006;
    tabentg[ 3][10][ 0] =  2.60858014891730109230E+006;
    tabentg[ 3][11][ 0] =  2.58198739656318211928E+006;
    tabentg[ 4][ 0][ 0] =  3.98784839662035694346E+006;
    tabentg[ 4][ 1][ 0] =  3.77635754210478533059E+006;
    tabentg[ 4][ 2][ 0] =  3.51789791043208679184E+006;
    tabentg[ 4][ 3][ 0] =  3.35684854448294593021E+006;
    tabentg[ 4][ 4][ 0] =  3.22237331442112335935E+006;
    tabentg[ 4][ 5][ 0] =  3.03767558785038720816E+006;
    tabentg[ 4][ 6][ 0] =  2.90726225009919889271E+006;
    tabentg[ 4][ 7][ 0] =  2.82811384822781337425E+006;
    tabentg[ 4][ 8][ 0] =  2.71134436910941870883E+006;
    tabentg[ 4][ 9][ 0] =  2.66559780960898287594E+006;
    tabentg[ 4][10][ 0] =  2.58615206360934022814E+006;
    tabentg[ 5][ 0][ 0] =  3.98509754539091885090E+006;
    tabentg[ 5][ 1][ 0] =  3.77278540319728432223E+006;
    tabentg[ 5][ 2][ 0] =  3.51289014005034742877E+006;
    tabentg[ 5][ 3][ 0] =  3.35048224774062633514E+006;
    tabentg[ 5][ 4][ 0] =  3.21424747987536853179E+006;
    tabentg[ 5][ 5][ 0] =  3.02464692848848318681E+006;
    tabentg[ 5][ 6][ 0] =  2.88604848782029002905E+006;
    tabentg[ 5][ 7][ 0] =  2.79748407805642439052E+006;
    tabentg[ 5][ 8][ 0] =  2.66642891078688064590E+006;
    tabentg[ 5][ 9][ 0] =  2.61736539279059320688E+006;
    tabentg[ 6][ 0][ 0] =  3.97640097055053291842E+006;
    tabentg[ 6][ 1][ 0] =  3.76148357130369357765E+006;
    tabentg[ 6][ 2][ 0] =  3.49698707169390469790E+006;
    tabentg[ 6][ 3][ 0] =  3.33007462136795744300E+006;
    tabentg[ 6][ 4][ 0] =  3.18763014066362706944E+006;
    tabentg[ 6][ 5][ 0] =  2.97853587750337226316E+006;
    tabentg[ 6][ 6][ 0] =  2.80298828443853789940E+006;
    tabentg[ 6][ 7][ 0] =  2.68574475614160299301E+006;
    tabentg[ 6][ 8][ 0] =  2.52546230677468702197E+006;
    tabentg[ 7][ 0][ 0] =  3.95717126664049457759E+006;
    tabentg[ 7][ 1][ 0] =  3.73641546847661491483E+006;
    tabentg[ 7][ 2][ 0] =  3.46118574507428379729E+006;
    tabentg[ 7][ 3][ 0] =  3.28259815515331458300E+006;
    tabentg[ 7][ 4][ 0] =  3.12171650581506639719E+006;
    tabentg[ 7][ 5][ 0] =  2.84244812027184758335E+006;
    tabentg[ 7][ 6][ 0] =  2.59417831099347816780E+006;
    tabentg[ 7][ 7][ 0] =  2.45213973343175742775E+006;
    tabentg[ 8][ 0][ 0] =  3.90869790000541741028E+006;
    tabentg[ 8][ 1][ 0] =  3.67240285060309385881E+006;
    tabentg[ 8][ 2][ 0] =  3.36457686016574036330E+006;
    tabentg[ 8][ 3][ 0] =  3.14240794951898232102E+006;
    tabentg[ 8][ 4][ 0] =  2.89474497843108884990E+006;
    tabentg[ 8][ 5][ 0] =  2.46091126350365346298E+006;
    tabentg[ 8][ 6][ 0] =  2.16390696532263653353E+006;
    tabentg[ 9][ 0][ 0] =  3.78484493337567104027E+006;
    tabentg[ 9][ 1][ 0] =  3.50190236055065412074E+006;
    tabentg[ 9][ 2][ 0] =  3.06708724457413749769E+006;
    tabentg[ 9][ 3][ 0] =  2.61169045174796506763E+006;
    tabentg[ 9][ 4][ 0] =  1.97800595516605814919E+006;
    tabentg[ 9][ 5][ 0] =  1.39376450072596827522E+006;
    tabentg[10][ 0][ 0] =  3.43374840990550303832E+006;
    tabentg[10][ 1][ 0] =  3.02490695981646468863E+006;
    tabentg[10][ 2][ 0] =  2.44776025615256559104E+006;
    tabentg[10][ 3][ 0] =  1.85372902841893257573E+006;
    tabentg[10][ 4][ 0] =  1.04210364191390445922E+006;
    tabentg[ 0][ 0][ 1] = -1.32477349858587931664E+001;
    tabentg[ 0][ 1][ 1] = -1.71986591769906915772E+001;
    tabentg[ 0][ 2][ 1] = -2.40836126954829552460E+001;
    tabentg[ 0][ 3][ 1] = -3.05257151512745927846E+001;
    tabentg[ 0][ 4][ 1] = -3.86743226784296112442E+001;
    tabentg[ 0][ 5][ 1] = -6.01109143223409319035E+001;
    tabentg[ 0][ 6][ 1] = -9.30674619121835604574E+001;
    tabentg[ 0][ 7][ 1] = -1.28435042398211209047E+002;
    tabentg[ 0][ 8][ 1] = -2.19148328110049931183E+002;
    tabentg[ 0][ 9][ 1] = -2.74306465159729782499E+002;
    tabentg[ 0][10][ 1] = -4.15278049110912093056E+002;
    tabentg[ 0][11][ 1] = -4.83341675729850692278E+002;
    tabentg[ 0][12][ 1] = -7.84053308160853134723E+002;
    tabentg[ 1][ 0][ 1] = -9.26681340044614643148E+001;
    tabentg[ 1][ 1][ 1] = -1.20305251628617938309E+002;
    tabentg[ 1][ 2][ 1] = -1.68467795136480418705E+002;
    tabentg[ 1][ 3][ 1] = -2.13538067634310181120E+002;
    tabentg[ 1][ 4][ 1] = -2.70562818663261623442E+002;
    tabentg[ 1][ 5][ 1] = -4.20682737967196771933E+002;
    tabentg[ 1][ 6][ 1] = -6.51729758325417606102E+002;
    tabentg[ 1][ 7][ 1] = -8.99933665379001354268E+002;
    tabentg[ 1][ 8][ 1] = -1.53740248317266355116E+003;
    tabentg[ 1][ 9][ 1] = -1.92548501858164695477E+003;
    tabentg[ 1][10][ 1] = -2.91859560066855283367E+003;
    tabentg[ 1][11][ 1] = -3.39864161867938810246E+003;
    tabentg[ 1][12][ 1] = -4.16839512446441767679E+003;
    tabentg[ 2][ 0][ 1] = -2.93040343284129619406E+002;
    tabentg[ 2][ 1][ 1] = -3.80438670474928926524E+002;
    tabentg[ 2][ 2][ 1] = -5.32758466583170161357E+002;
    tabentg[ 2][ 3][ 1] = -6.75343802773471793444E+002;
    tabentg[ 2][ 4][ 1] = -8.55872574319707041468E+002;
    tabentg[ 2][ 5][ 1] = -1.33196198943026683992E+003;
    tabentg[ 2][ 6][ 1] = -2.06673564373658973636E+003;
    tabentg[ 2][ 7][ 1] = -2.85810813735056535734E+003;
    tabentg[ 2][ 8][ 1] = -4.89748240065721256542E+003;
    tabentg[ 2][ 9][ 1] = -6.14278661223703329597E+003;
    tabentg[ 2][10][ 1] = -9.33960980642778849870E+003;
    tabentg[ 2][11][ 1] = -9.08377536797233005927E+003;
    tabentg[ 2][12][ 1] = -1.20568220410522717430E+004;
    tabentg[ 3][ 0][ 1] = -9.26654918432080762614E+002;
    tabentg[ 3][ 1][ 1] = -1.20305422525451012916E+003;
    tabentg[ 3][ 2][ 1] = -1.68490118362007387987E+003;
    tabentg[ 3][ 3][ 1] = -2.13641770665099602411E+003;
    tabentg[ 3][ 4][ 1] = -2.70935569245737633537E+003;
    tabentg[ 3][ 5][ 1] = -4.22888481038803547563E+003;
    tabentg[ 3][ 6][ 1] = -6.59475183864538666967E+003;
    tabentg[ 3][ 7][ 1] = -9.16352016215087132878E+003;
    tabentg[ 3][ 8][ 1] = -1.58526812667321701156E+004;
    tabentg[ 3][ 9][ 1] = -1.99750983244549352094E+004;
    tabentg[ 3][10][ 1] = -2.30509481959100230597E+004;
    tabentg[ 3][11][ 1] = -2.56230046261913303169E+004;
    tabentg[ 4][ 0][ 1] = -2.93014803211265234495E+003;
    tabentg[ 4][ 1][ 1] = -3.80443526870755567870E+003;
    tabentg[ 4][ 2][ 1] = -5.33001645610118248442E+003;
    tabentg[ 4][ 3][ 1] = -6.76455345375775686989E+003;
    tabentg[ 4][ 4][ 1] = -8.59841516176516233827E+003;
    tabentg[ 4][ 5][ 1] = -1.35531032151177532796E+004;
    tabentg[ 4][ 6][ 1] = -2.14841151486415219551E+004;
    tabentg[ 4][ 7][ 1] = -3.03078432435899230768E+004;
    tabentg[ 4][ 8][ 1] = -5.39854438948846873245E+004;
    tabentg[ 4][ 9][ 1] = -5.39734423381038141088E+004;
    tabentg[ 4][10][ 1] = -6.66613930359334335662E+004;
    tabentg[ 5][ 0][ 1] = -9.26427112284631584771E+003;
    tabentg[ 5][ 1][ 1] = -1.20320115754938469763E+004;
    tabentg[ 5][ 2][ 1] = -1.68795073913381456805E+004;
    tabentg[ 5][ 3][ 1] = -2.14980246947063787957E+004;
    tabentg[ 5][ 4][ 1] = -2.75619596531814022455E+004;
    tabentg[ 5][ 5][ 1] = -4.49772125666264910251E+004;
    tabentg[ 5][ 6][ 1] = -7.51569173049178789370E+004;
    tabentg[ 5][ 7][ 1] = -1.10820315119216073072E+005;
    tabentg[ 5][ 8][ 1] = -1.25676389395267557120E+005;
    tabentg[ 5][ 9][ 1] = -1.38956224935454869410E+005;
    tabentg[ 6][ 0][ 1] = -2.92871422110774619796E+004;
    tabentg[ 6][ 1][ 1] = -3.80889901819162623724E+004;
    tabentg[ 6][ 2][ 1] = -5.37877984649906502455E+004;
    tabentg[ 6][ 3][ 1] = -6.96194604326253029285E+004;
    tabentg[ 6][ 4][ 1] = -9.24937268908403930254E+004;
    tabentg[ 6][ 5][ 1] = -1.69464932503145246301E+005;
    tabentg[ 6][ 6][ 1] = -3.25087069932358223014E+005;
    tabentg[ 6][ 7][ 1] = -3.36136972540069487877E+005;
    tabentg[ 6][ 8][ 1] = -4.38190026653506909497E+005;
    tabentg[ 7][ 0][ 1] = -7.36099297957523958758E+004;
    tabentg[ 7][ 1][ 1] = -9.62291845838674780680E+004;
    tabentg[ 7][ 2][ 1] = -1.39235479188555793371E+005;
    tabentg[ 7][ 3][ 1] = -1.89414720231925253756E+005;
    tabentg[ 7][ 4][ 1] = -2.73781438116773322690E+005;
    tabentg[ 7][ 5][ 1] = -6.37509843864739639685E+005;
    tabentg[ 7][ 6][ 1] = -7.18962797292940318584E+005;
    tabentg[ 7][ 7][ 1] = -8.31888141009158222005E+005;
    tabentg[ 8][ 0][ 1] = -1.86247880095038184663E+005;
    tabentg[ 8][ 1][ 1] = -2.48803623452488565817E+005;
    tabentg[ 8][ 2][ 1] = -3.92093720179372467101E+005;
    tabentg[ 8][ 3][ 1] = -6.02255580573899904266E+005;
    tabentg[ 8][ 4][ 1] = -1.12650760046413401142E+006;
    tabentg[ 8][ 5][ 1] = -1.27017443997623119503E+006;
    tabentg[ 8][ 6][ 1] = -1.43239393106126808561E+006;
    tabentg[ 9][ 0][ 1] = -4.78412323265573999379E+005;
    tabentg[ 9][ 1][ 1] = -6.74562869461943395436E+005;
    tabentg[ 9][ 2][ 1] = -1.28458059961046115495E+006;
    tabentg[ 9][ 3][ 1] = -3.16724768985062744468E+006;
    tabentg[ 9][ 4][ 1] = -3.45718751586101902649E+006;
    tabentg[ 9][ 5][ 1] = -4.06555937391219427809E+006;
    tabentg[10][ 0][ 1] = -1.05094385816645156592E+006;
    tabentg[10][ 1][ 1] = -1.11904387122455728240E+006;
    tabentg[10][ 2][ 1] = -7.56915826508124475367E+005;
    tabentg[10][ 3][ 1] = -1.62884528838215395808E+005;
    tabentg[10][ 4][ 1] = -6.54734741630255710334E+005;
    tabentg[ 0][ 0][ 2] =  1.47817924328558612615E+006;
    tabentg[ 0][ 1][ 2] =  1.43744729096597549506E+006;
    tabentg[ 0][ 2][ 2] =  1.38459163152905274183E+006;
    tabentg[ 0][ 3][ 2] =  1.35092707118482864462E+006;
    tabentg[ 0][ 4][ 2] =  1.32302607989458041266E+006;
    tabentg[ 0][ 5][ 2] =  1.28608640671582054347E+006;
    tabentg[ 0][ 6][ 2] =  1.26188945295818126760E+006;
    tabentg[ 0][ 7][ 2] =  1.24839150499627296813E+006;
    tabentg[ 0][ 8][ 2] =  1.23071625503042363562E+006;
    tabentg[ 0][ 9][ 2] =  1.22452599260627035983E+006;
    tabentg[ 0][10][ 2] =  1.21456553053938923404E+006;
    tabentg[ 0][11][ 2] =  1.21134798112467466854E+006;
    tabentg[ 0][12][ 2] =  1.20251683699850365520E+006;
    tabentg[ 1][ 0][ 2] =  1.47824151150185544975E+006;
    tabentg[ 1][ 1][ 2] =  1.43752811833979981020E+006;
    tabentg[ 1][ 2][ 2] =  1.38471165235648024827E+006;
    tabentg[ 1][ 3][ 2] =  1.35109907623011199757E+006;
    tabentg[ 1][ 4][ 2] =  1.32328939495986746624E+006;
    tabentg[ 1][ 5][ 2] =  1.28668166768826683983E+006;
    tabentg[ 1][ 6][ 2] =  1.26309863486600527540E+006;
    tabentg[ 1][ 7][ 2] =  1.25030363567222375423E+006;
    tabentg[ 1][ 8][ 2] =  1.23448442281157267280E+006;
    tabentg[ 1][ 9][ 2] =  1.22943056999708036892E+006;
    tabentg[ 1][10][ 2] =  1.22236283454779442400E+006;
    tabentg[ 1][11][ 2] =  1.22053099732198519632E+006;
    tabentg[ 1][12][ 2] =  1.21824411677364422940E+006;
    tabentg[ 2][ 0][ 2] =  1.47839861305063753389E+006;
    tabentg[ 2][ 1][ 2] =  1.43773204763180133887E+006;
    tabentg[ 2][ 2][ 2] =  1.38501450607694406062E+006;
    tabentg[ 2][ 3][ 2] =  1.35153324256056058221E+006;
    tabentg[ 2][ 4][ 2] =  1.32395445975431357510E+006;
    tabentg[ 2][ 5][ 2] =  1.28818778703314508311E+006;
    tabentg[ 2][ 6][ 2] =  1.26616432515557785518E+006;
    tabentg[ 2][ 7][ 2] =  1.25515939116605441086E+006;
    tabentg[ 2][ 8][ 2] =  1.24408215276671946049E+006;
    tabentg[ 2][ 9][ 2] =  1.24194046404788712971E+006;
    tabentg[ 2][10][ 2] =  1.24230958321907022037E+006;
    tabentg[ 2][11][ 2] =  1.23518535546307149343E+006;
    tabentg[ 2][12][ 2] =  1.23203048947499692440E+006;
    tabentg[ 3][ 0][ 2] =  1.47889542716224794276E+006;
    tabentg[ 3][ 1][ 2] =  1.43837698091087280773E+006;
    tabentg[ 3][ 2][ 2] =  1.38597266088048112579E+006;
    tabentg[ 3][ 3][ 2] =  1.35290816759711643681E+006;
    tabentg[ 3][ 4][ 2] =  1.32606460633918200620E+006;
    tabentg[ 3][ 5][ 2] =  1.29299178899074438959E+006;
    tabentg[ 3][ 6][ 2] =  1.27600237003213050775E+006;
    tabentg[ 3][ 7][ 2] =  1.27081662317252578214E+006;
    tabentg[ 3][ 8][ 2] =  1.27530194426178699359E+006;
    tabentg[ 3][ 9][ 2] =  1.28280013906019879505E+006;
    tabentg[ 3][10][ 2] =  1.27668446627273247577E+006;
    tabentg[ 3][11][ 2] =  1.25595861507192929275E+006;
    tabentg[ 4][ 0][ 2] =  1.48046665460041468032E+006;
    tabentg[ 4][ 1][ 2] =  1.44041699884398793802E+006;
    tabentg[ 4][ 2][ 2] =  1.38900736136394948699E+006;
    tabentg[ 4][ 3][ 2] =  1.35727687101880018599E+006;
    tabentg[ 4][ 4][ 2] =  1.33281143795478902757E+006;
    tabentg[ 4][ 5][ 2] =  1.30861575632238062099E+006;
    tabentg[ 4][ 6][ 2] =  1.30861288027562643401E+006;
    tabentg[ 4][ 7][ 2] =  1.32347844651226070710E+006;
    tabentg[ 4][ 8][ 2] =  1.38305108180222497322E+006;
    tabentg[ 4][ 9][ 2] =  1.35626385834363452159E+006;
    tabentg[ 4][10][ 2] =  1.33681227723747235723E+006;
    tabentg[ 5][ 0][ 2] =  1.48543723741029365920E+006;
    tabentg[ 5][ 1][ 2] =  1.44687487229654355906E+006;
    tabentg[ 5][ 2][ 2] =  1.39865995105441473424E+006;
    tabentg[ 5][ 3][ 2] =  1.37133335938244778663E+006;
    tabentg[ 5][ 4][ 2] =  1.35499431454988219775E+006;
    tabentg[ 5][ 5][ 2] =  1.36288003878097725101E+006;
    tabentg[ 5][ 6][ 2] =  1.42826247845554677770E+006;
    tabentg[ 5][ 7][ 2] =  1.52423174577393126674E+006;
    tabentg[ 5][ 8][ 2] =  1.50244648216010397300E+006;
    tabentg[ 5][ 9][ 2] =  1.43548872719842661172E+006;
    tabentg[ 6][ 0][ 2] =  1.50118384972776775248E+006;
    tabentg[ 6][ 1][ 2] =  1.46739989412676193751E+006;
    tabentg[ 6][ 2][ 2] =  1.42999411145756649785E+006;
    tabentg[ 6][ 3][ 2] =  1.41910830374435521662E+006;
    tabentg[ 6][ 4][ 2] =  1.43624257683417177759E+006;
    tabentg[ 6][ 5][ 2] =  1.59261121365371323191E+006;
    tabentg[ 6][ 6][ 2] =  1.99681430494577880017E+006;
    tabentg[ 6][ 7][ 2] =  1.89831952551717637107E+006;
    tabentg[ 6][ 8][ 2] =  1.80335367002360778861E+006;
    tabentg[ 7][ 0][ 2] =  1.53626554882081551477E+006;
    tabentg[ 7][ 1][ 2] =  1.51377649932139297016E+006;
    tabentg[ 7][ 2][ 2] =  1.50637955500786099583E+006;
    tabentg[ 7][ 3][ 2] =  1.55083103929845243692E+006;
    tabentg[ 7][ 4][ 2] =  1.69273410709662479348E+006;
    tabentg[ 7][ 5][ 2] =  2.52921712992360256612E+006;
    tabentg[ 7][ 6][ 2] =  2.44612974940243782476E+006;
    tabentg[ 7][ 7][ 2] =  2.27275987058828631416E+006;
    tabentg[ 8][ 0][ 2] =  1.62734325122457579710E+006;
    tabentg[ 8][ 1][ 2] =  1.64128515457327873446E+006;
    tabentg[ 8][ 2][ 2] =  1.76449373410690738820E+006;
    tabentg[ 8][ 3][ 2] =  2.08425505279741412960E+006;
    tabentg[ 8][ 4][ 2] =  3.14638101610660273582E+006;
    tabentg[ 8][ 5][ 2] =  3.05124348285675933585E+006;
    tabentg[ 8][ 6][ 2] =  2.90074646065059397370E+006;
    tabentg[ 9][ 0][ 2] =  1.88229906022444926202E+006;
    tabentg[ 9][ 1][ 2] =  2.05583894770053401589E+006;
    tabentg[ 9][ 2][ 2] =  2.94178719774645287544E+006;
    tabentg[ 9][ 3][ 2] =  7.50218492351750936359E+006;
    tabentg[ 9][ 4][ 2] =  5.17150500812062807381E+006;
    tabentg[ 9][ 5][ 2] =  3.17480148388065490872E+006;
    tabentg[10][ 0][ 2] =  2.60113015676527610049E+006;
    tabentg[10][ 1][ 2] =  3.13383917064896645024E+006;
    tabentg[10][ 2][ 2] =  3.60680052731625083834E+006;
    tabentg[10][ 3][ 2] =  6.52160915271159354597E+006;
    tabentg[10][ 4][ 2] =  9.71089857738896831870E+006;
    // !;
    // ! liquid table [density];
    // !;
    tabdenl[ 0][ 0][ 0] =  9.99841676114285178301E+002;
    tabdenl[ 0][ 1][ 0] =  9.95370961989185161656E+002;
    tabdenl[ 0][ 2][ 0] =  9.83005720176408203770E+002;
    tabdenl[ 1][ 0][ 0] =  9.99844667110630439311E+002;
    tabdenl[ 1][ 1][ 0] =  9.95376995022998698914E+002;
    tabdenl[ 1][ 2][ 0] =  9.83024785551920103899E+002;
    tabdenl[ 1][ 3][ 0] =  9.73808766901105855140E+002;
    tabdenl[ 2][ 0][ 0] =  9.99852213141648689998E+002;
    tabdenl[ 2][ 1][ 0] =  9.95383701200031396183E+002;
    tabdenl[ 2][ 2][ 0] =  9.83036069643778319005E+002;
    tabdenl[ 2][ 3][ 0] =  9.73845402063280630500E+002;
    tabdenl[ 2][ 4][ 0] =  9.41089225602210376564E+002;
    tabdenl[ 3][ 0][ 0] =  9.99876074654955573351E+002;
    tabdenl[ 3][ 1][ 0] =  9.95404906413740263815E+002;
    tabdenl[ 3][ 2][ 0] =  9.83057261169481989782E+002;
    tabdenl[ 3][ 3][ 0] =  9.73882037225455405860E+002;
    tabdenl[ 3][ 4][ 0] =  9.41221811177524045888E+002;
    tabdenl[ 3][ 5][ 0] =  9.12452440052769361500E+002;
    tabdenl[ 4][ 0][ 0] =  9.99951520065697650352E+002;
    tabdenl[ 4][ 1][ 0] =  9.95471947396969426336E+002;
    tabdenl[ 4][ 2][ 0] =  9.83124255293924761645E+002;
    tabdenl[ 4][ 3][ 0] =  9.73950826367300919628E+002;
    tabdenl[ 4][ 4][ 0] =  9.41354396752837715212E+002;
    tabdenl[ 4][ 5][ 0] =  9.12691236521278710825E+002;
    tabdenl[ 4][ 6][ 0] =  8.58284683136555145211E+002;
    tabdenl[ 5][ 0][ 0] =  1.00018998621438083774E+003;
    tabdenl[ 5][ 1][ 0] =  9.95683792034678958771E+002;
    tabdenl[ 5][ 2][ 0] =  9.83335916112309519121E+002;
    tabdenl[ 5][ 3][ 0] =  9.74168140336080341513E+002;
    tabdenl[ 5][ 4][ 0] =  9.41603168502390872163E+002;
    tabdenl[ 5][ 5][ 0] =  9.13051330867454453255E+002;
    tabdenl[ 5][ 6][ 0] =  8.61747637346084275123E+002;
    tabdenl[ 5][ 7][ 0] =  7.97878206851450045178E+002;
    tabdenl[ 6][ 0][ 0] =  1.00094295067599193771E+003;
    tabdenl[ 6][ 1][ 0] =  9.96352138538192093620E+002;
    tabdenl[ 6][ 2][ 0] =  9.84003327970179043405E+002;
    tabdenl[ 6][ 3][ 0] =  9.74853199149043462057E+002;
    tabdenl[ 6][ 4][ 0] =  9.42386702840913130785E+002;
    tabdenl[ 6][ 5][ 0] =  9.13952052220719679099E+002;
    tabdenl[ 6][ 6][ 0] =  8.64447012677121961133E+002;
    tabdenl[ 6][ 7][ 0] =  8.02918411454867737120E+002;
    tabdenl[ 6][ 8][ 0] =  6.95007412907598450147E+002;
    tabdenl[ 7][ 0][ 0] =  1.00260172155163706975E+003;
    tabdenl[ 7][ 1][ 0] =  9.97821614491530795021E+002;
    tabdenl[ 7][ 2][ 0] =  9.85468871462965807950E+002;
    tabdenl[ 7][ 3][ 0] =  9.76356548728260463577E+002;
    tabdenl[ 7][ 4][ 0] =  9.44102524595787713224E+002;
    tabdenl[ 7][ 5][ 0] =  9.15920046759583442508E+002;
    tabdenl[ 7][ 6][ 0] =  8.47352168353659521927E+002;
    tabdenl[ 7][ 7][ 0] =  8.08933253636413496679E+002;
    tabdenl[ 7][ 8][ 0] =  6.98275061010562694719E+002;
    tabdenl[ 7][ 9][ 0] =  4.69460801747971459008E+002;
    tabdenl[ 8][ 0][ 0] =  1.00673111307154715632E+003;
    tabdenl[ 8][ 1][ 0] =  1.00146434608260199184E+003;
    tabdenl[ 8][ 2][ 0] =  9.89090981399345650971E+002;
    tabdenl[ 8][ 3][ 0] =  9.80066639790722433645E+002;
    tabdenl[ 8][ 4][ 0] =  9.48316332492425203782E+002;
    tabdenl[ 8][ 5][ 0] =  9.20728316018421310218E+002;
    tabdenl[ 8][ 6][ 0] =  8.54313437459811098051E+002;
    tabdenl[ 8][ 7][ 0] =  7.99007195899779617321E+002;
    tabdenl[ 8][ 8][ 0] =  7.02529917047158960486E+002;
    tabdenl[ 8][ 9][ 0] =  4.87347675017029075661E+002;
    tabdenl[ 8][10][ 0] = -1.68237148055064153596E+002;
    tabdenl[ 9][ 0][ 0] =  1.01686114577972273310E+003;
    tabdenl[ 9][ 1][ 0] =  1.01033682622326148248E+003;
    tabdenl[ 9][ 2][ 0] =  9.97852531467949120270E+002;
    tabdenl[ 9][ 3][ 0] =  9.89010396159166930374E+002;
    tabdenl[ 9][ 4][ 0] =  9.58362744205881085691E+002;
    tabdenl[ 9][ 5][ 0] =  9.32062960811865764299E+002;
    tabdenl[ 9][ 6][ 0] =  8.70070289225560372870E+002;
    tabdenl[ 9][ 7][ 0] =  8.20534533634368358435E+002;
    tabdenl[ 9][ 8][ 0] =  7.42150226285701251072E+002;
    tabdenl[ 9][ 9][ 0] =  5.11342274388889791226E+002;
    tabdenl[ 9][10][ 0] = -1.15960107872597831147E+002;
    tabdenl[10][ 0][ 0] =  1.04537484401641495424E+003;
    tabdenl[10][ 1][ 0] =  1.03567632027373429082E+003;
    tabdenl[10][ 2][ 0] =  1.02251668001294433452E+003;
    tabdenl[10][ 3][ 0] =  1.01398598936062421672E+003;
    tabdenl[10][ 4][ 0] =  9.85732973612952605436E+002;
    tabdenl[10][ 5][ 0] =  9.62217361368654792386E+002;
    tabdenl[10][ 6][ 0] =  9.08866394382883413527E+002;
    tabdenl[10][ 7][ 0] =  8.68804886085472617197E+002;
    tabdenl[10][ 8][ 0] =  8.11754667604265364389E+002;
    tabdenl[10][ 9][ 0] =  7.00609416630898294898E+002;
    tabdenl[10][10][ 0] =  5.02861515601031101141E+002;
    tabdenl[ 0][ 0][ 1] =  1.14879276830054280305E-003;
    tabdenl[ 0][ 1][ 1] =  7.14145994053965067760E-003;
    tabdenl[ 0][ 2][ 1] =  2.25681837162683418541E-002;
    tabdenl[ 1][ 0][ 1] =  8.03577468547832414469E-003;
    tabdenl[ 1][ 1][ 1] =  7.14145994064169665344E-003;
    tabdenl[ 1][ 2][ 1] =  2.25681837163297788207E-002;
    tabdenl[ 1][ 3][ 1] =  7.32703243495375911021E-002;
    tabdenl[ 2][ 0][ 1] =  2.54107715351941272630E-002;
    tabdenl[ 2][ 1][ 1] =  2.25824705716299971303E-002;
    tabdenl[ 2][ 2][ 1] =  2.25681837165306459214E-002;
    tabdenl[ 2][ 3][ 1] =  7.32703243495638478766E-002;
    tabdenl[ 2][ 4][ 1] =  2.65171150627321106086E-001;
    tabdenl[ 3][ 0][ 1] =  8.03501228795839711472E-002;
    tabdenl[ 3][ 1][ 1] =  7.14039580204523238738E-002;
    tabdenl[ 3][ 2][ 1] =  7.13569492265149346588E-002;
    tabdenl[ 3][ 3][ 1] =  7.32703243495375911021E-002;
    tabdenl[ 3][ 4][ 1] =  2.65171150627356189133E-001;
    tabdenl[ 3][ 5][ 1] =  4.77592937018684882133E-001;
    tabdenl[ 4][ 0][ 1] =  2.54031478087475615091E-001;
    tabdenl[ 4][ 1][ 1] =  2.25718375849041308356E-001;
    tabdenl[ 4][ 2][ 1] =  2.25551444830033143152E-001;
    tabdenl[ 4][ 3][ 1] =  2.31590173272320815068E-001;
    tabdenl[ 4][ 4][ 1] =  2.65171150627321106086E-001;
    tabdenl[ 4][ 5][ 1] =  4.77592937018712415664E-001;
    tabdenl[ 4][ 6][ 1] =  7.54170928663467776687E+000;
    tabdenl[ 5][ 0][ 1] =  8.02738923746154964256E-001;
    tabdenl[ 5][ 1][ 1] =  7.12978929545705297777E-001;
    tabdenl[ 5][ 2][ 1] =  7.12268972177674375423E-001;
    tabdenl[ 5][ 3][ 1] =  7.31246625681438944255E-001;
    tabdenl[ 5][ 4][ 1] =  8.36921381543873077469E-001;
    tabdenl[ 5][ 5][ 1] =  9.62784447684257305156E-001;
    tabdenl[ 5][ 6][ 1] =  6.31010755148184188101E+000;
    tabdenl[ 5][ 7][ 1] =  6.28050953697895941730E+000;
    tabdenl[ 6][ 0][ 1] =  2.53269276931364739625E+000;
    tabdenl[ 6][ 1][ 1] =  2.24666019074349021523E+000;
    tabdenl[ 6][ 2][ 1] =  2.24261592722780855880E+000;
    tabdenl[ 6][ 3][ 1] =  2.30145830493868830402E+000;
    tabdenl[ 6][ 4][ 1] =  2.63053228171153952175E+000;
    tabdenl[ 6][ 5][ 1] =  3.02179809268060983740E+000;
    tabdenl[ 6][ 6][ 1] =  4.48739377266890215878E+000;
    tabdenl[ 6][ 7][ 1] =  1.38803088766918083508E+001;
    tabdenl[ 6][ 8][ 1] =  5.70110042333055755392E+000;
    tabdenl[ 7][ 0][ 1] =  6.32965171703985785001E+000;
    tabdenl[ 7][ 1][ 1] =  5.60004511489133083302E+000;
    tabdenl[ 7][ 2][ 1] =  5.58018146924313640511E+000;
    tabdenl[ 7][ 3][ 1] =  5.72167927757470273775E+000;
    tabdenl[ 7][ 4][ 1] =  6.52096512082009471811E+000;
    tabdenl[ 7][ 5][ 1] =  7.46792326412315432549E+000;
    tabdenl[ 7][ 6][ 1] =  1.09596412029170320324E+001;
    tabdenl[ 7][ 7][ 1] =  1.61939020310369876654E+001;
    tabdenl[ 7][ 8][ 1] =  1.06371400914906644175E+001;
    tabdenl[ 7][ 9][ 1] =  2.94478679156363014613E+001;
    tabdenl[ 8][ 0][ 1] =  1.56946530978901446218E+001;
    tabdenl[ 8][ 1][ 1] =  1.38076769582130491898E+001;
    tabdenl[ 8][ 2][ 1] =  1.37007825708102650708E+001;
    tabdenl[ 8][ 3][ 1] =  1.40191363226526402030E+001;
    tabdenl[ 8][ 4][ 1] =  1.58686263867755901202E+001;
    tabdenl[ 8][ 5][ 1] =  1.80432711615992715792E+001;
    tabdenl[ 8][ 6][ 1] =  2.57811028772457255798E+001;
    tabdenl[ 8][ 7][ 1] =  3.65758079791612047416E+001;
    tabdenl[ 8][ 8][ 1] =  7.66170941573065817920E+001;
    tabdenl[ 8][ 9][ 1] =  5.99864984296517889106E+001;
    tabdenl[ 8][10][ 1] =  1.30692600456165791911E+002;
    tabdenl[ 9][ 0][ 1] =  3.80604600062484124123E+001;
    tabdenl[ 9][ 1][ 1] =  3.32437791298052118805E+001;
    tabdenl[ 9][ 2][ 1] =  3.26777670271780138478E+001;
    tabdenl[ 9][ 3][ 1] =  3.32801132680568088063E+001;
    tabdenl[ 9][ 4][ 1] =  3.71094504177208150963E+001;
    tabdenl[ 9][ 5][ 1] =  4.15602751651046702364E+001;
    tabdenl[ 9][ 6][ 1] =  5.63147331139905631403E+001;
    tabdenl[ 9][ 7][ 1] =  7.43543483911833220645E+001;
    tabdenl[ 9][ 8][ 1] =  1.22879860360924638485E+002;
    tabdenl[ 9][ 9][ 1] =  7.92356032292025929564E+002;
    tabdenl[ 9][10][ 1] =  2.01685805240938589122E+003;
    tabdenl[10][ 0][ 1] =  9.22223105797211815116E+001;
    tabdenl[10][ 1][ 1] =  8.55079917149064527848E+001;
    tabdenl[10][ 2][ 1] =  8.26608811108622347774E+001;
    tabdenl[10][ 3][ 1] =  8.32792646964506388940E+001;
    tabdenl[10][ 4][ 1] =  8.98691795152487173937E+001;
    tabdenl[10][ 5][ 1] =  9.76135247781182329163E+001;
    tabdenl[10][ 6][ 1] =  1.20206244263811754536E+002;
    tabdenl[10][ 7][ 1] =  1.42506179943847001823E+002;
    tabdenl[10][ 8][ 1] =  1.85131759158555269096E+002;
    tabdenl[10][ 9][ 1] =  3.13347859110689341833E+002;
    tabdenl[10][10][ 1] =  7.01958186153046312938E+002;
    tabdenl[ 0][ 0][ 2] =  4.45924771514057667332E+001;
    tabdenl[ 0][ 1][ 2] = -2.36057864346099734121E+002;
    tabdenl[ 0][ 2][ 2] = -3.13508438443986790389E+002;
    tabdenl[ 1][ 0][ 2] =  4.45793678930188690401E+001;
    tabdenl[ 1][ 1][ 2] = -2.03089243949629349117E+002;
    tabdenl[ 1][ 2][ 2] = -3.45897843653863276359E+002;
    tabdenl[ 1][ 3][ 2] = -3.91383648411277590640E+002;
    tabdenl[ 2][ 0][ 2] =  4.45462925541691134868E+001;
    tabdenl[ 2][ 1][ 2] = -2.03095789064314089956E+002;
    tabdenl[ 2][ 2][ 2] = -3.37350719403072901059E+002;
    tabdenl[ 2][ 3][ 2] = -3.97902687036743145654E+002;
    tabdenl[ 2][ 4][ 2] = -5.35321713563548428283E+002;
    tabdenl[ 3][ 0][ 2] =  4.44416872763449504191E+001;
    tabdenl[ 3][ 1][ 2] = -2.03116494516760042188E+002;
    tabdenl[ 3][ 2][ 2] = -3.37334757364669769686E+002;
    tabdenl[ 3][ 3][ 2] = -3.95258831264797379390E+002;
    tabdenl[ 3][ 4][ 2] = -5.35231939331537091675E+002;
    tabdenl[ 3][ 5][ 2] = -6.20164491381503466982E+002;
    tabdenl[ 4][ 0][ 2] =  4.41107781622374872654E+001;
    tabdenl[ 4][ 1][ 2] = -2.03182050650180542561E+002;
    tabdenl[ 4][ 2][ 2] = -3.37284366229860097519E+002;
    tabdenl[ 4][ 3][ 2] = -3.95167059579771830613E+002;
    tabdenl[ 4][ 4][ 2] = -5.29287159756982191539E+002;
    tabdenl[ 4][ 5][ 2] = -6.21843773237356344907E+002;
    tabdenl[ 4][ 6][ 2] = -5.01098505085111298740E+002;
    tabdenl[ 5][ 0][ 2] =  4.30631871358023303742E+001;
    tabdenl[ 5][ 1][ 2] = -2.03390152214599822855E+002;
    tabdenl[ 5][ 2][ 2] = -3.37125865819048954108E+002;
    tabdenl[ 5][ 3][ 2] = -3.94877871534876533133E+002;
    tabdenl[ 5][ 4][ 2] = -5.28676293893361616938E+002;
    tabdenl[ 5][ 5][ 2] = -6.18798596189174986648E+002;
    tabdenl[ 5][ 6][ 2] = -4.40101166893800666458E+002;
    tabdenl[ 5][ 7][ 2] = -1.58428648745276700538E+003;
    tabdenl[ 6][ 0][ 2] =  3.97393166587990478433E+001;
    tabdenl[ 6][ 1][ 2] = -2.04056094621344186635E+002;
    tabdenl[ 6][ 2][ 2] = -3.36633075587151893160E+002;
    tabdenl[ 6][ 3][ 2] = -3.93973478753038136801E+002;
    tabdenl[ 6][ 4][ 2] = -5.26764512839330222960E+002;
    tabdenl[ 6][ 5][ 2] = -6.15941556019060271865E+002;
    tabdenl[ 6][ 6][ 2] = -8.35767785035532824622E+002;
    tabdenl[ 6][ 7][ 2] = -1.11442559760326935248E+003;
    tabdenl[ 6][ 8][ 2] = -1.58334936607846316292E+003;
    tabdenl[ 7][ 0][ 2] =  3.23374953792573194278E+001;
    tabdenl[ 7][ 1][ 2] = -2.05569906563652523346E+002;
    tabdenl[ 7][ 2][ 2] = -3.35588013713117106818E+002;
    tabdenl[ 7][ 3][ 2] = -3.92027016944038848578E+002;
    tabdenl[ 7][ 4][ 2] = -5.22641849988111061975E+002;
    tabdenl[ 7][ 5][ 2] = -6.09807606369129757695E+002;
    tabdenl[ 7][ 6][ 2] = -8.20807238342756136262E+002;
    tabdenl[ 7][ 7][ 2] = -1.05332763749405467024E+003;
    tabdenl[ 7][ 8][ 2] = -1.71312717815221594719E+003;
    tabdenl[ 7][ 9][ 2] = -1.80709219511841843087E+003;
    tabdenl[ 8][ 0][ 2] =  1.35412477430552886659E+001;
    tabdenl[ 8][ 1][ 2] = -2.09607166713159017490E+002;
    tabdenl[ 8][ 2][ 2] = -3.33220031507122996572E+002;
    tabdenl[ 8][ 3][ 2] = -3.87443421097564112188E+002;
    tabdenl[ 8][ 4][ 2] = -5.12877025502507422061E+002;
    tabdenl[ 8][ 5][ 2] = -5.95422098721931092769E+002;
    tabdenl[ 8][ 6][ 2] = -7.87150635868297172237E+002;
    tabdenl[ 8][ 7][ 2] = -9.81542412168711962295E+002;
    tabdenl[ 8][ 8][ 2] = -1.53664566195446059282E+003;
    tabdenl[ 8][ 9][ 2] = -1.77385036927830674358E+003;
    tabdenl[ 8][10][ 2] = -5.22652147345627145114E+003;
    tabdenl[ 9][ 0][ 2] = -3.30719888348416333201E+001;
    tabdenl[ 9][ 1][ 2] = -2.20972694652353936817E+002;
    tabdenl[ 9][ 2][ 2] = -3.28706802299384094113E+002;
    tabdenl[ 9][ 3][ 2] = -3.77618914280568560571E+002;
    tabdenl[ 9][ 4][ 2] = -4.91504694560997279495E+002;
    tabdenl[ 9][ 5][ 2] = -5.64588873889331239297E+002;
    tabdenl[ 9][ 6][ 2] = -7.20940575380045174825E+002;
    tabdenl[ 9][ 7][ 2] = -8.56146645095937628867E+002;
    tabdenl[ 9][ 8][ 2] = -1.13422137260273325410E+003;
    tabdenl[ 9][ 9][ 2] = -3.34918516957548126811E+003;
    tabdenl[ 9][10][ 2] = -3.34918516957548126811E+003;
    tabdenl[10][ 0][ 2] = -1.16379028529749604104E+002;
    tabdenl[10][ 1][ 2] = -2.60206000276574684449E+002;
    tabdenl[10][ 2][ 2] = -3.24034998459412918237E+002;
    tabdenl[10][ 3][ 2] = -3.58170659481389918710E+002;
    tabdenl[10][ 4][ 2] = -4.44480470004698361208E+002;
    tabdenl[10][ 5][ 2] = -4.99354207681284549381E+002;
    tabdenl[10][ 6][ 2] = -6.01532180791107862206E+002;
    tabdenl[10][ 7][ 2] = -6.68431715864277407491E+002;
    tabdenl[10][ 8][ 2] = -7.61626969045263308544E+002;
    tabdenl[10][ 9][ 2] = -9.21386621565653399557E+002;
    tabdenl[10][10][ 2] = -1.09595621887325864918E+003;
    // !;
    // ! liquid table [enthalpy];
    // !;
    tabentl[ 0][ 0][ 0] =  3.55070305991915165578E+003;
    tabentl[ 0][ 1][ 0] =  1.30147664337938505923E+005;
    tabentl[ 0][ 2][ 0] =  2.51867219722832407570E+005;
    tabentl[ 1][ 0][ 0] =  3.55681105853479266443E+003;
    tabentl[ 1][ 1][ 0] =  1.30160023263288065209E+005;
    tabentl[ 1][ 2][ 0] =  2.51903254383869643789E+005;
    tabentl[ 1][ 3][ 0] =  3.19611921986047353130E+005;
    tabentl[ 2][ 0][ 0] =  3.57222103754644467699E+003;
    tabentl[ 2][ 1][ 0] =  1.30173758878638691385E+005;
    tabentl[ 2][ 2][ 0] =  2.51924581969312042929E+005;
    tabentl[ 2][ 3][ 0] =  3.19676491029577504378E+005;
    tabentl[ 2][ 4][ 0] =  5.11069518243931815960E+005;
    tabentl[ 3][ 0][ 0] =  3.62094949652990771938E+003;
    tabentl[ 3][ 1][ 0] =  1.30217193854756784276E+005;
    tabentl[ 3][ 2][ 0] =  2.51964630670714483131E+005;
    tabentl[ 3][ 3][ 0] =  3.19741060073107655626E+005;
    tabentl[ 3][ 4][ 0] =  5.11247232516376592685E+005;
    tabentl[ 3][ 5][ 0] =  6.48613705863796873018E+005;
    tabentl[ 4][ 0][ 0] =  3.77502067971407495861E+003;
    tabentl[ 4][ 1][ 0] =  1.30354538773781285272E+005;
    tabentl[ 4][ 2][ 0] =  2.52091273997555486858E+005;
    tabentl[ 4][ 3][ 0] =  3.19862308383457420859E+005;
    tabentl[ 4][ 4][ 0] =  5.11424946788821369410E+005;
    tabentl[ 4][ 5][ 0] =  6.49192005162840010598E+005;
    tabentl[ 4][ 6][ 0] =  8.85371834611636819318E+005;
    tabentl[ 5][ 0][ 0] =  4.26201919496796017484E+003;
    tabentl[ 5][ 1][ 0] =  1.30788776152640071814E+005;
    tabentl[ 5][ 2][ 0] =  2.52491737414172064746E+005;
    tabentl[ 5][ 3][ 0] =  3.20245740562132734340E+005;
    tabentl[ 5][ 4][ 0] =  5.11758737358004902489E+005;
    tabentl[ 5][ 5][ 0] =  6.49724899841658654623E+005;
    tabentl[ 5][ 6][ 0] =  8.80089316162579809316E+005;
    tabentl[ 5][ 7][ 0] =  1.06996759110597078688E+006;
    tabentl[ 6][ 0][ 0] =  5.79987024429833581962E+003;
    tabentl[ 6][ 1][ 0] =  1.32161100351978064282E+005;
    tabentl[ 6][ 2][ 0] =  2.53757932035736448597E+005;
    tabentl[ 6][ 3][ 0] =  3.21458369403604068793E+005;
    tabentl[ 6][ 4][ 0] =  5.12815312588947068434E+005;
    tabentl[ 6][ 5][ 0] =  6.50641933271307847463E+005;
    tabentl[ 6][ 6][ 0] =  8.78709606850157841109E+005;
    tabentl[ 6][ 7][ 0] =  1.07314309475553780794E+006;
    tabentl[ 6][ 8][ 0] =  1.30658672469775914215E+006;
    tabentl[ 7][ 0][ 0] =  9.18848711018657013483E+003;
    tabentl[ 7][ 1][ 0] =  1.35190816143560601631E+005;
    tabentl[ 7][ 2][ 0] =  2.56556597485491773114E+005;
    tabentl[ 7][ 3][ 0] =  3.24140152304698189255E+005;
    tabentl[ 7][ 4][ 0] =  5.15156949170183681417E+005;
    tabentl[ 7][ 5][ 0] =  6.52680847517268033698E+005;
    tabentl[ 7][ 6][ 0] =  9.29733102896619704552E+005;
    tabentl[ 7][ 7][ 0] =  1.07289193811953766271E+006;
    tabentl[ 7][ 8][ 0] =  1.36925940723576489836E+006;
    tabentl[ 7][ 9][ 0] =  1.70949521732148854062E+006;
    tabentl[ 8][ 0][ 0] =  1.76293371415661349602E+004;
    tabentl[ 8][ 1][ 0] =  1.42772837891154485987E+005;
    tabentl[ 8][ 2][ 0] =  2.63579915925686538685E+005;
    tabentl[ 8][ 3][ 0] =  3.30879086232746602036E+005;
    tabentl[ 8][ 4][ 0] =  5.21069889371886791196E+005;
    tabentl[ 8][ 5][ 0] =  6.57866367080797092058E+005;
    tabentl[ 8][ 6][ 0] =  9.32547840609947685152E+005;
    tabentl[ 8][ 7][ 0] =  1.12236972535291942768E+006;
    tabentl[ 8][ 8][ 0] =  1.39230887141490494832E+006;
    tabentl[ 8][ 9][ 0] =  1.90185452714877459221E+006;
    tabentl[ 8][10][ 0] =  3.12198868620389886200E+006;
    tabentl[ 9][ 0][ 0] =  3.83939298022216753452E+004;
    tabentl[ 9][ 1][ 0] =  1.61637583874034025939E+005;
    tabentl[ 9][ 2][ 0] =  2.81173690532040840480E+005;
    tabentl[ 9][ 3][ 0] =  3.47813281751838920172E+005;
    tabentl[ 9][ 4][ 0] =  5.36090265962692210451E+005;
    tabentl[ 9][ 5][ 0] =  6.71237164655239554122E+005;
    tabentl[ 9][ 6][ 0] =  9.40855529253966524266E+005;
    tabentl[ 9][ 7][ 0] =  1.12430748846596828662E+006;
    tabentl[ 9][ 8][ 0] =  1.37409055166748445481E+006;
    tabentl[ 9][ 9][ 0] =  1.92041335258496832103E+006;
    tabentl[ 9][10][ 0] =  3.15698578387281624600E+006;
    tabentl[10][ 0][ 0] =  9.91218365225149027538E+004;
    tabentl[10][ 1][ 0] =  2.18089867653321212856E+005;
    tabentl[10][ 2][ 0] =  3.34776154583011346404E+005;
    tabentl[10][ 3][ 0] =  3.99812529782366298605E+005;
    tabentl[10][ 4][ 0] =  5.83322352422941708937E+005;
    tabentl[10][ 5][ 0] =  7.14487439424493582919E+005;
    tabentl[10][ 6][ 0] =  9.73173345922558219172E+005;
    tabentl[10][ 7][ 0] =  1.14529814202432567254E+006;
    tabentl[10][ 8][ 0] =  1.36983512384330132045E+006;
    tabentl[10][ 9][ 0] =  1.76463494048038101755E+006;
    tabentl[10][10][ 0] =  2.39776025615256559104E+006;
    tabentl[ 0][ 0][ 1] =  2.34640080632802572325E+000;
    tabentl[ 0][ 1][ 1] =  1.46295832280238826684E+001;
    tabentl[ 0][ 2][ 1] =  4.26551708848049528910E+001;
    tabentl[ 1][ 0][ 1] =  1.64130197512821638384E+001;
    tabentl[ 1][ 1][ 1] =  1.46295832280247406487E+001;
    tabentl[ 1][ 2][ 1] =  4.26551708847876369646E+001;
    tabentl[ 1][ 3][ 1] =  1.29138087060324977529E+002;
    tabentl[ 2][ 0][ 1] =  5.19014134077147133439E+001;
    tabentl[ 2][ 1][ 1] =  4.62623674558354665010E+001;
    tabentl[ 2][ 2][ 1] =  4.26551708848089248249E+001;
    tabentl[ 2][ 3][ 1] =  1.29138087060280014384E+002;
    tabentl[ 2][ 4][ 1] =  3.55428544889553904795E+002;
    tabentl[ 3][ 0][ 1] =  1.64115557229019060514E+002;
    tabentl[ 3][ 1][ 1] =  1.46290083230143977744E+002;
    tabentl[ 3][ 2][ 1] =  1.34886580501789438813E+002;
    tabentl[ 3][ 3][ 1] =  1.29138087060324977529E+002;
    tabentl[ 3][ 4][ 1] =  3.55428544889552995301E+002;
    tabentl[ 3][ 5][ 1] =  1.15659859808622968558E+003;
    tabentl[ 4][ 0][ 1] =  5.18867729801333894102E+002;
    tabentl[ 4][ 1][ 1] =  4.62566175560716658310E+002;
    tabentl[ 4][ 2][ 1] =  4.26539667887753637388E+002;
    tabentl[ 4][ 3][ 1] =  4.08376349019367864912E+002;
    tabentl[ 4][ 4][ 1] =  3.55428544889553904795E+002;
    tabentl[ 4][ 5][ 1] =  1.15659859808632063505E+003;
    tabentl[ 4][ 6][ 1] = -1.40980466131589164434E+004;
    tabentl[ 5][ 0][ 1] =  1.63969149611182660919E+003;
    tabentl[ 5][ 1][ 1] =  1.46232557533684712325E+003;
    tabentl[ 5][ 2][ 1] =  1.34874478825353821776E+003;
    tabentl[ 5][ 3][ 1] =  1.29145717068278850093E+003;
    tabentl[ 5][ 4][ 1] =  1.12449899465162661727E+003;
    tabentl[ 5][ 5][ 1] =  9.74980117188255462679E+002;
    tabentl[ 5][ 6][ 1] = -7.03202718306912265689E+003;
    tabentl[ 5][ 7][ 1] =  1.33299061882684472948E+004;
    tabentl[ 6][ 0][ 1] =  5.17403668639109309879E+003;
    tabentl[ 6][ 1][ 1] =  4.61990091277418014215E+003;
    tabentl[ 6][ 2][ 1] =  4.26416750962912828982E+003;
    tabentl[ 6][ 3][ 1] =  4.08449732722253793327E+003;
    tabentl[ 6][ 4][ 1] =  3.56125144630133809187E+003;
    tabentl[ 6][ 5][ 1] =  3.09407375159983166668E+003;
    tabentl[ 6][ 6][ 1] =  1.51318993338125005721E+003;
    tabentl[ 6][ 7][ 1] = -6.27891590000363066792E+002;
    tabentl[ 6][ 8][ 1] =  1.53809521271466044709E+005;
    tabentl[ 7][ 0][ 1] =  1.29348575527188058913E+004;
    tabentl[ 7][ 1][ 1] =  1.15802514885778819007E+004;
    tabentl[ 7][ 2][ 1] =  1.07056673297028173693E+004;
    tabentl[ 7][ 3][ 1] =  1.02625056600810275995E+004;
    tabentl[ 7][ 4][ 1] =  8.97371066067694846424E+003;
    tabentl[ 7][ 5][ 1] =  7.83035272023947345588E+003;
    tabentl[ 7][ 6][ 1] =  4.02847946638463781710E+003;
    tabentl[ 7][ 7][ 1] = -1.53742956828077808495E+003;
    tabentl[ 7][ 8][ 1] =  1.59553891418562707258E+005;
    tabentl[ 7][ 9][ 1] =  4.80898274568215187173E+005;
    tabentl[ 8][ 0][ 1] =  3.21035488778323015140E+004;
    tabentl[ 8][ 1][ 1] =  2.89329713832534544053E+004;
    tabentl[ 8][ 2][ 1] =  2.68538698407316987868E+004;
    tabentl[ 8][ 3][ 1] =  2.57903890192376566119E+004;
    tabentl[ 8][ 4][ 1] =  2.27045752001975088206E+004;
    tabentl[ 8][ 5][ 1] =  2.00062565000800641428E+004;
    tabentl[ 8][ 6][ 1] =  1.13830595377034333069E+004;
    tabentl[ 8][ 7][ 1] = -2.94861325100947738065E+002;
    tabentl[ 8][ 8][ 1] = -4.43065705228624501615E+004;
    tabentl[ 8][ 9][ 1] =  9.99938217965988791548E+005;
    tabentl[ 8][10][ 1] =  2.75375577104191249236E+006;
    tabentl[ 9][ 0][ 1] =  7.82922152716022683308E+004;
    tabentl[ 9][ 1][ 1] =  7.16847060483829845907E+004;
    tabentl[ 9][ 2][ 1] =  6.71732742588774417527E+004;
    tabentl[ 9][ 3][ 1] =  6.47920861417693959083E+004;
    tabentl[ 9][ 4][ 1] =  5.78710700259620571160E+004;
    tabentl[ 9][ 5][ 1] =  5.19827684369029229856E+004;
    tabentl[ 9][ 6][ 1] =  3.45907583780517379637E+004;
    tabentl[ 9][ 7][ 1] =  1.42070689170968817052E+004;
    tabentl[ 9][ 8][ 1] = -4.08389592472731455928E+004;
    tabentl[ 9][ 9][ 1] = -9.07144090785020147450E+005;
    tabentl[ 9][10][ 1] = -2.57877028269732603803E+006;
    tabentl[10][ 0][ 1] =  2.09883708732048864476E+005;
    tabentl[10][ 1][ 1] =  1.96073088502896396676E+005;
    tabentl[10][ 2][ 1] =  1.88436234532422735356E+005;
    tabentl[10][ 3][ 1] =  1.83767445541751774726E+005;
    tabentl[10][ 4][ 1] =  1.69318794585564231966E+005;
    tabentl[10][ 5][ 1] =  1.57335800421289488440E+005;
    tabentl[10][ 6][ 1] =  1.26431671964990280685E+005;
    tabentl[10][ 7][ 1] =  9.78265914475614990806E+004;
    tabentl[10][ 8][ 1] =  4.41301656026673445012E+004;
    tabentl[10][ 9][ 1] = -1.16532576427351508755E+005;
    tabentl[10][10][ 1] = -7.56915826508124475367E+005;
    tabentl[ 0][ 0][ 2] =  2.72755715707732131705E+006;
    tabentl[ 0][ 1][ 2] =  2.69415747966997744516E+006;
    tabentl[ 0][ 2][ 2] =  2.71560053743641357869E+006;
    tabentl[ 1][ 0][ 2] =  2.72753672440574830398E+006;
    tabentl[ 1][ 1][ 2] =  2.70431250133059499785E+006;
    tabentl[ 1][ 2][ 2] =  2.70649777069524815306E+006;
    tabentl[ 1][ 3][ 2] =  2.71019563747897557914E+006;
    tabentl[ 2][ 0][ 2] =  2.72748517796706361696E+006;
    tabentl[ 2][ 1][ 2] =  2.70428577389534702525E+006;
    tabentl[ 2][ 2][ 2] =  2.70842042127879476175E+006;
    tabentl[ 2][ 3][ 2] =  2.71173230354244960472E+006;
    tabentl[ 2][ 4][ 2] =  2.74106049458729987964E+006;
    tabentl[ 3][ 0][ 2] =  2.72732221420891815796E+006;
    tabentl[ 3][ 1][ 2] =  2.70420127357883332297E+006;
    tabentl[ 3][ 2][ 2] =  2.70835082135777734220E+006;
    tabentl[ 3][ 3][ 2] =  2.71419534947413485497E+006;
    tabentl[ 3][ 4][ 2] =  2.74182095945090288296E+006;
    tabentl[ 3][ 5][ 2] =  2.77490487779489625245E+006;
    tabentl[ 4][ 0][ 2] =  2.72680728145253797993E+006;
    tabentl[ 4][ 1][ 2] =  2.70393425264530861750E+006;
    tabentl[ 4][ 2][ 2] =  2.70813084966762922704E+006;
    tabentl[ 4][ 3][ 2] =  2.71398208907442353666E+006;
    tabentl[ 4][ 4][ 2] =  2.74908374678218178451E+006;
    tabentl[ 4][ 5][ 2] =  2.78372984253584081307E+006;
    tabentl[ 4][ 6][ 2] =  2.09098283958586771041E+006;
    tabentl[ 5][ 0][ 2] =  2.72518295878006564453E+006;
    tabentl[ 5][ 1][ 2] =  2.70309178299907129258E+006;
    tabentl[ 5][ 2][ 2] =  2.70743646258377050981E+006;
    tabentl[ 5][ 3][ 2] =  2.71330884783840458840E+006;
    tabentl[ 5][ 4][ 2] =  2.74830063214136054739E+006;
    tabentl[ 5][ 5][ 2] =  2.79622494968071533367E+006;
    tabentl[ 5][ 6][ 2] =  1.95845856571497349069E+006;
    tabentl[ 5][ 7][ 2] =  4.05987027559694182128E+006;
    tabentl[ 6][ 0][ 2] =  2.72008672960285143927E+006;
    tabentl[ 6][ 1][ 2] =  2.70044688724222825840E+006;
    tabentl[ 6][ 2][ 2] =  2.70525281948862923309E+006;
    tabentl[ 6][ 3][ 2] =  2.71119128408289561048E+006;
    tabentl[ 6][ 4][ 2] =  2.74584103675205027685E+006;
    tabentl[ 6][ 5][ 2] =  2.79299241891608014703E+006;
    tabentl[ 6][ 6][ 2] =  2.97180664453957695514E+006;
    tabentl[ 6][ 7][ 2] =  3.19090295626485766843E+006;
    tabentl[ 6][ 8][ 2] =  2.64518779229067638516E+006;
    tabentl[ 7][ 0][ 2] =  2.70903518326858244836E+006;
    tabentl[ 7][ 1][ 2] =  2.69470220547721302137E+006;
    tabentl[ 7][ 2][ 2] =  2.70048982494070706889E+006;
    tabentl[ 7][ 3][ 2] =  2.70656989610573602840E+006;
    tabentl[ 7][ 4][ 2] =  2.74049121196681587026E+006;
    tabentl[ 7][ 5][ 2] =  2.78599912624898133799E+006;
    tabentl[ 7][ 6][ 2] =  2.95592288839830597863E+006;
    tabentl[ 7][ 7][ 2] =  3.19494220909927412868E+006;
    tabentl[ 7][ 8][ 2] =  4.21424451880640815943E+006;
    tabentl[ 7][ 9][ 2] =  1.02015255943549377844E+006;
    tabentl[ 8][ 0][ 2] =  2.68257042627881933004E+006;
    tabentl[ 8][ 1][ 2] =  2.68089580387374106795E+006;
    tabentl[ 8][ 2][ 2] =  2.68891392162158293650E+006;
    tabentl[ 8][ 3][ 2] =  2.69531898534297756851E+006;
    tabentl[ 8][ 4][ 2] =  2.72756008478184556589E+006;
    tabentl[ 8][ 5][ 2] =  2.76929938568564504385E+006;
    tabentl[ 8][ 6][ 2] =  2.91968363703364692628E+006;
    tabentl[ 8][ 7][ 2] =  3.11794583502138452604E+006;
    tabentl[ 8][ 8][ 2] =  3.76516234704183042049E+006;
    tabentl[ 8][ 9][ 2] =  4.07400158732539461926E+006;
    tabentl[ 8][10][ 2] =  8.95466001497171446681E+006;
    tabentl[ 9][ 0][ 2] =  2.62361736679531121626E+006;
    tabentl[ 9][ 1][ 2] =  2.65005351663256855682E+006;
    tabentl[ 9][ 2][ 2] =  2.66217180431233346462E+006;
    tabentl[ 9][ 3][ 2] =  2.66915655714380973950E+006;
    tabentl[ 9][ 4][ 2] =  2.69786884236024599522E+006;
    tabentl[ 9][ 5][ 2] =  2.73193206949582090601E+006;
    tabentl[ 9][ 6][ 2] =  2.84590598988474858925E+006;
    tabentl[ 9][ 7][ 2] =  2.98009006403027614579E+006;
    tabentl[ 9][ 8][ 2] =  3.30921597725148545578E+006;
    tabentl[ 9][ 9][ 2] =  6.60209520175039023161E+006;
    tabentl[ 9][10][ 2] =  6.60209520175038930029E+006;
    tabentl[10][ 0][ 2] =  2.46839779717351822183E+006;
    tabentl[10][ 1][ 2] =  2.58508409011563798413E+006;
    tabentl[10][ 2][ 2] =  2.59847174474308127537E+006;
    tabentl[10][ 3][ 2] =  2.60456252955194935203E+006;
    tabentl[10][ 4][ 2] =  2.62457866954716714099E+006;
    tabentl[10][ 5][ 2] =  2.64399090196688612923E+006;
    tabentl[10][ 6][ 2] =  2.70025880520274629816E+006;
    tabentl[10][ 7][ 2] =  2.75904964575967518613E+006;
    tabentl[10][ 8][ 2] =  2.85460753088593855500E+006;
    tabentl[10][ 9][ 2] =  3.28403848968593776226E+006;
    tabentl[10][10][ 2] =  3.60680052731625083834E+006;

    return;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
/* 
int main ()
{
    double a = 300;
    double psatwater,viscwater,hfgwater; 

    psatwater = psatwaterf(a);
    viscwater = viscwaterf(100,200);
    // eptab()

    hfgwater =hfgwaterf(100);

    printf("DONE\n");
    return 0;
} 
*/

/* 
int main ()
{
    wrsttabset();
    printf("wrsttabset.DONE\n");
    double a, b, d,e,g, h,i,j;
    int c,f; 
    a = 500.0;
    b = 100;
    c = 1;

    d = 10;
    e = 20;
    g = 51;

    h = 25;
    i = 22;
    j= 33;

    f = 1;


    f = wrsteamtab( a,b,c, d,e,g,  h,i,j);


    printf("DONE1\n");
    printf("DONE2\n");
    return 0;
} 
*/

#endif