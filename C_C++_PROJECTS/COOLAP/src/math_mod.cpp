#ifndef MATH_MOD_H
#define MATH_MOD_H

#include <vector>
#include <iostream>
// #include <string>
// #include <cmath>
#include <math.h>

using namespace std;



// void gauss (double **a[][3], int n, int ist=0,int neps=0,int iscl=0, double eps=0.0)
void gauss (vector<vector<long double>> &a, int n, int &ist,int &neps,int iscl) 
{
    int i, j, k, max;
    int l; // whang
    double t, xx, xmax, xeps=1e-20;



    // ! scaling if wanted
    // if(present(iscl)) {
        if (iscl == 1) {
            for (i=1; i<=n; i++) // do i=1, n
            {
                max=1;
                xmax=fabs(a[i-1][max-1]);
                for (j=2; j<=n; j++) // do j=2, n
                {
                    xx=fabs(a[i-1][j-1]);
                    if( xx > xmax ) {
                        xmax = xx; max = j; 
                    }
                }
                if( xmax < xeps ) {  // ! nearly 0 coeffs.
                    ist = 1;
                    neps = i;
                    return;
                }
                // a(i,1:n+1) = a(i,1:n+1)/xmax;
                for (l=0;l<=n;l++)
                {
                    a[i-1][l] = a[i-1][l]/xmax;
                }
            }
        }
    // end if

    // !
    // ! pivot and gaussian elimination on the matrix
    // !

    for (k=1; k<=n-1; k++) // do k=1, n-1
    {
        max=k;
        xmax = fabs(a[k-1][k-1]);
        for (i=k+1; i<=n; i++) // do i=k+1, n
        {
            xx = fabs(a[i-1][k-1]);
            if( xx > xmax ) {
                xmax = xx; max=i;
            }
        }
        if( xmax < xeps ) {  // ! nearly 0 pivot => singlar
            ist = 1;
            neps = k;
            return;
        }
        if( max != k ) {
            for (j=k;j<=n+1;j++) // do j=k, n+1
            {
                t = a[k-1][j-1];
                a[k-1][j-1] = a[max-1][j-1];
                a[max-1][j-1] = t;
            }
        }
        for (i=k+1; i<=n; i++) // do i=k+1, n // ! do the elimination
        {
            t=a[i-1][k-1]/a[k-1][k-1];
            for (j=k+1;j<=n+1; j++) { // do j=k+1, n+1  // ! j<=k elemects are already zeroed out // !
                a[i-1][j-1]=a[i-1][j-1]-a[k-1][j-1]*t;
            }
        }
    }
    // !
    // ! solve it
    // !
    for (k=n;k>0;k--) // do k=n, 1, -1
    {
        t = a[k-1][n];
        for (j=k+1; j<=n; j++) // do j=k+1, n
        {
            t = t-a[k-1][j-1]*a[j-1][n];
        }
        a[k-1][n] = t/a[k-1][k-1];
    }
    ist = 0;


}




#endif