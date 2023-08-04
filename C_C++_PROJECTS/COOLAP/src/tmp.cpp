// #ifndef 
// #define 



#include <iostream>
#include <fstream>

#include <string>

// namespace fs = files


using namespace std;

#include <math.h>
#include <random>
#include <vector>

// #include <filesystem>
// namespace fs = std::filesystem;



// #include <direct.h>
// #include <sys/stat.h>

// #include <filesystem>
// namespace fs = std::filesystem;



/*###############################################################################*/ 

int main()
{

    // int i ;
    // const int nthexmax = 3;

    // double l;
    // l = log(-2);

    // vector<double> thex = vector<double>(10000000,0.0);
    // vector<double> thex1 = vector<double>(10000000,0.0);
    // vector<double> thex2 = vector<double>(10000000,0.0);
    // vector<double> thex3 = vector<double>(10000000,0.0);
    // static double ffpsz[100000];
    // static double ffpsz1[100000];
    // static double ffpsz2[100000];
    // static double ffpsz3[100000];
    // static int mat [10000][1000];
    // vector<double> t, x;
    // string tetete="test";

    // #include "opr1k.dat"

    // vector<double> thex  = vector<double>(nthexmax, 0.0);
    // vector<double> qhex  = vector<double>(nthexmax, 0.0);
    // double qhex2[nthexmax]={};

    // for (i=1; i<=nthexmax;i++)
    // {
    //     thex[i]=222;
    //     qhex[i-1]=111.0;
    //     qhex2[i]=111;
    // }

    // if (thex[nthexmax]>100||qhex[nthexmax]>100){
    //     printf("her44e\n");
    // }
    // qhex2[5]=24;


    system("mkdir -p mbrq");
    
    ofstream oj1, oj2;
    oj1.open("mbrq/oj1", ios::out);
    oj1 << "TEST" << endl;
    oj1.close();

    oj2.open("mbrq/oj2", ios::out);
    oj2 << "TEST2" << endl;
    oj2.close();



    printf("################\n");
    printf("DONE\n");

/* 
    ifstream ifs;
    string str, mode;

    // 267 ! (before Do-loop) open    :: ("i",ombrq,md)
    // 267 ! (before Do-loop) headers :: ("w",ombrq,md, time, istep)
    // 606 ! (end of Do-loop) write   :: ("w",ombrq,md, time, istep)
    // 625 ! (out of Do-loop) close   :: ("c",ombrq,md)

    mode ="i";

    // outdbrcool(mode,odir, md, time,step,ist)
    // ! can be called only by the "mode" as i (initialize), w (write) or c (close)

    ofstream fout;
    fout.open("zTest_ujh", ios::out);
    // fout << "HEADER~~~~" << endl;
    fout.close();
    fout.open("zTest_ujh", ios::app);
    

    mode="w";

    for (int i = 0; i<=10; i++)
    {
        if (mode=="i"){
            printf ("Init \n");
            // fout.open("zTest_ujh");
            // fout.open("zTest_ujh", ios::app);
            // fout << "HEADER~~~~" << endl;
            // header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

        }
        else if (mode=="w"){
            // ofstream fout("zTest_ujh");

            printf ("Write \n");
            // 1688 ! write hex data when a section of time is passed
            // fout.open("zTest_ujh",ios::app);
            // 1700 ! history data
            fout << i << "\t" << i << endl;

        }
        // else {
        //     printf ("Close \n");
        // }
    }
    fout.close();
 */

/* 
    double lj=3000.0;

    // header
    // contents - tab separated...

    ofstream fout("ztest2.txt");
    fout << lj << "\t" << "Tab!?" << endl;
    fout.close();
*/


/* 
    ifs.open("ztest.txt");
    if(!ifs.is_open())
    {
        printf("failed\n");
        exit(1);
    }
    else
    {
        printf("write\n");
        getline(ifs, str);
        cout << str << endl;
        ifs.close();
    }
 */

    return 0;
}

// #endif

// #include "dbrcool.cpp"
// #include "meltprop.cpp"
// #include "consts.cpp"
// #include "coolpropdef.cpp"
// #include "wrsttab.cpp"
// #include "htpar.cpp"



//     // vector<double> c;
//     int n;
//     n=4;
//     vector<vector<double>> b = vector<vector<double>>(n,vector<double>(n+1,0.0));
//     b[0][0]=1;
//     b[2][1]=44;


//     vector<double> a1 = vector<double>(3,0.0);


//     vectorTest(b, n);

//     fill(a1.begin(),a1.end(),442);
//     fill(b.begin(),b.end(),vector<double>(n+1,4422));

//     printf("hello \n");

    


// // void vectorTest (const int n, double *a[][n])
// void vectorTest (vector<vector<double>> a, int n)
// {
//     for (int i=0; i<n; i++)
//     {
//         cout << a[i][0] << ";";
//         printf("test\n");
//     }
// }