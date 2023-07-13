#include <iostream>

int main()
{
    for (int kn = 0; kn < numnod; kn++)
    {
        if (srctot <= 0)
        {
            if (nactiv <= 0 || htot <= 0)
            {
                recon = 0;
            }
            else
            {
                recon = denmlt / visrex;
            }
        }
        else
        {
            for (int ktp = 0; ktp < 16; ktp++)
            {
                void index(ktp, kj);
                rvle = romliq(kj);
                if (ktp == 16)
                {
                    rvle = rslagl * roc;
                    dsrdt(kn) = dsrdt(kn) + (double source(ktp, kn) * double area(kn)) / rvle;
                }
                else
                {
                    dsrdt(kn) = dsrdt(kn) + (double source(ktp, kn) * double area(kn)) / rvle;
                }
            }
        }
    }
    dtvel = dtime / nvelp;
    for (int kitv = 0; kitv < nvelp; kitv++)
    {
        if (kitv != 1)
        {
            for (int ksv = 1; ksv < numnod; ksv++)
            {
                velold(ksv) = vel(ksv);
                htmold(ksv) = htmp(ksv);
            }
            velold(numnod + 1) = vel(numnod + 1);
        }
        else
            ;
        niter = 0;

        while (1)
        {
            niter++;
            if(niter > nitmax){
                std::cout << "warning t = " << time << '\n' << std::endl;
                break;
            }
            else{
                nup = numnod + 1;
                
            }






        }
    }

    jsend = 1;
    if (ngeom != 1)
    {
        for (int i = 0; i < numnod; i++)
        {
            if (iflga(i) == 1)
                jsend = 1;
            else
                ;
        }
    }
    else
        ;

    for (int kn = 0; kn < numnod; kn++)
    {
        if ((nactiv(kn) != 0) && (htot(kn) <= 0) && (htmp(kn) > 0) && (srctot(kn) <= 0))
        {
            void IO34();
        }
        else if ((nactiv(kn) == 0) || (htot(kn) > 0) || (htmp(kn) <= 0))
        {
            continue;
        }
        else
        {
            void IO33();
        }
    }
    system("pause");
    return 0;
}