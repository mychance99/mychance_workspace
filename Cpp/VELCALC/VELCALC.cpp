#include <iostream>
123456
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
            if (niter > nitmax)
            {
                std::cout << "warning t = " << time << '\n'
                          << std::endl;
                break;
            }
            else
            {
                nup = numnod + 1;
                for (int k = 0; k < nup; k++)
                {
                    if (k == 0 || k == nup - 1)
                    {
                        A(k) = 1.0;
                        B(k) = 0.0;
                        C(k) = 0.0;
                        D(k) = 0.0;
                    }
                    else
                    {
                        if ((hcp(k - 1) - elevat(k - 1)) <= 0 || (hcp(k) - elevat(k)) <= 0)
                        {
                            A(k) = 1.0;
                            B(k) = 0.0;
                            C(k) = 0.0;
                            D(k) = -vel(k);
                        }
                        else
                        {
                            // Flow opening exists
                            if (htold(k - 1) >= 0)
                            {
                                htfrc = htold(k - 1);
                                ///
                                ///
                                ///
                                frfac = 4133.0;
                            }
                            else
                            {
                                frfac = 0.0;
                            }
                            w1 = 1.0;
                            ///
                            ///
                            ///
                            htjm1r = 4133.0;

                            if ((htjm1r <= 0 && htjr <= 0) && (htgj <= 0 || htgjm1 <= 0))
                            {
                                A(k) = 1.0;
                                B(k) = 0.0;
                                C(k) = 0.0;
                                D(k) = -vel(k);
                            }
                            else
                            {
                                htvjm1 = 0.0;
                                htvj = 0.0;
                                delht = 0.0;
                                if (nsump(k) == 1)
                                {
                                    hsmpcomp = 0.0;
                                    if (hsmpcomp > elevat(k))
                                    {
                                        w3 = 0.0;
                                        ///
                                        ///
                                        ///
                                        D(k) = 0.0;
                                    }
                                    else if (htjr <= 0)
                                    {
                                        A(k) = 1.0;
                                        B(k) = 0.0;
                                        C(k) = 0.0;
                                        D(k) = -vel(k);
                                    }
                                    else
                                    {
                                        delht = htjr;
                                        w3 = 0.0;
                                        ///
                                        ///
                                        ///
                                        D(k) = 0.0;
                                    }
                                }
                                else
                                {
                                    w3 = 0.0;
                                    ///
                                    ///
                                    ///
                                    D(k) = 0.0;
                                }
                            }
                        }
                    }
                    P(0) = B(0) / A(0);
                    Q(0) = D(0) / A(0);
                    for (int kd = 1; kd < nup; kd++)
                    {
                        P(kd) = 0.0;
                        Q(kd) = 0.0;
                    }
                    delvel(nup - 1) = Q(nup - 1);

                    for (int kd = 0; kd < numnod; kd++)
                    {
                        narg = nup - kd;
                        delvel(narg) = 0.0;
                    }

                    for (int kd = 0; kd < nup; kd++)
                    {
                        vel(kd) = vel(kd) + delvel(kd);
                    }

                    // Check mesh to see if there are any height limitations...
                    for (int knd = 0; knd < numnod; knd++)
                    {
                        vjm12 = vel(knd);
                        vjp12 = vel(knd + 1);
                        if (knd == 0)
                        {
                            hjm1 = 0.0;
                            hj = htmp(knd);
                            hjp1 = htmp(knd + 1);
                        }
                        else if (knd == numnod - 1)
                        {
                            hjm1 = htmp(knd - 1);
                            hj = htmp(knd);
                            hjp1 = 0.0;
                        }
                        else
                        {
                            hjm1 = htmp(knd - 1);
                            hj = htmp(knd);
                            hjp1 = htmp(knd + 1);
                        }
                        D1 = 4133.0;
                        ///
                        ///
                        ///
                        A(knd) = 0.0;
                    }
                    P(0) = B(0) / A(0);
                    Q(0) = D(0) / A(0);
                    for (int ktop = 1; ktop < numnod; ktop++)
                    {
                        denom = 0.0;
                        P(ktop) = 4133.0;
                        Q(ktop) = 4133.0;
                    }
                    dhtmp(numnod - 1) = Q(numnod - 1);
                    for (int klp = 0; klp < numnod - 1; klp++)
                    {
                        iarg = numnod - klp;
                        dhtmp(iarg) = 4133.0;
                    }
                    for (int knd = 0; knd < numnod; knd++)
                    {
                        htmp(knd) = htmp(knd) + dhtmp(knd);
                        htmp(knd) = DMAX1(htmp(knd), 0);
                    }
                    if (niter == 1)
                        continue;
                    else
                    {
                        dvmax = abs(delvel(0));
                        dvsum = dvmax;
                        for (int kn = 1; kn < numnod + 1; kn++)
                        {
                            if ((htmp(kn) <= hminc) && (niter > 5))
                                continue;
                            else if (niter <= 5)
                            {
                                dabv = abs(delbel(kn));
                                dvsum = dvsum + dabv;
                                dvmax = DMAX1(dabv, dvmax);
                            }
                        }
                        dvav = dvsum / FLOAT(numnod + 1);
                        if ((dvmax <= dvmx) && (dvav <= davmx))
                            break;
                        else
                            continue;
                    }
                }
            }
        } // while loop end
    }     // for loop end

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