#ifndef MISCMOD_H
#define MISCMOD_H

using namespace std;

#include <vector>


/*###############################################################################*/
void swap (double &x1, double &x2)
{
    double x;
    x = x1;
    x1 = x2;
    x2 = x;
    return;
}
/*###############################################################################*/

/*###############################################################################*/ 
void split(vector<double> &x, int istarr, int low, int high, int &mid, vector<int> &order)
{
    int left, right;
    double pivot;
    left = low;
    right = high;
    pivot = x[low-1];

    while (left >= right) {
        while (left >= right || x[right-1] < pivot){
            right = right-1;
        }
        while (left >= right || x[left] > pivot) {
            left = left+1;
        }
         if(left < right) {
            swap(x[left-1], x[right-1]);
            swap(order[left-1], order[right-1]);
         }
    }
    swap(x[low-1], x[right-1]);
    swap(order[low-1], order[right-1]);
    mid = right;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void qsort(vector<double> &x, int istarr, int first, int last, vector<int> &order)
{
    int mid;
    if (first < last) {
        split(x, istarr, first, last, mid, order);
        qsort(x, istarr, first, mid-1, order);
        qsort(x, istarr, mid+1, last, order);
    }
    return;
}

/*###############################################################################*/ 




#endif