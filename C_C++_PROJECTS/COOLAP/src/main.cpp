#include <iostream>
#include <math.h>
#include <random>
#include <fstream>
#include <algorithm>

using namespace std;

#include "hotcan.cpp"
#include "hotbody.cpp"
#include "inttable_mod.cpp"
#include "dbrcool.cpp"
#include "dbrcinteg.cpp"

// #include "consts.cpp"
// #include "wrsttab.cpp"
// #include "meltprop.cpp"
// #include "coolpropdef.cpp"
// #include "htpar.cpp"


    // as GLOBAL VARIABLES... some array size is too big... bd, bdbk... 
    /* STATEMENT 26-65 */
    const int  ncmax=5, nbmax = 20, nmvmax=10, uoc0 = 30, uob0 = 40;


    /*##################################################################*/ 
    // TYPE = struct
    // hotcan, htbody, etc.
    meio xio;
    meltdebr md, mdbk;
    htcan cn, cnbk;
    htbody bd[nbmax], bdbk[nbmax]; 
    timetab tdt, tdtout;



//  int main(int argc, char* argv[])
int main()
{

    /*##################################################################*/ 
    int  uerr=0, uin=10, udb=11, /* uoc[ncmax],  uob[ncmax][nbmax], */ uob1, narg, 
        ist=0, idbg=0, icmma=0, itu=0,
        ncan=1, nbody[ncmax]={}, nmv=0, nbtot, i,j, ioutcount, ngood, is,is1,is2,is3,
        ostep=500, istep = 0;
    // double free or corruption (out) >> string-related?
    // string strtu="s";
    // string ombrq="mbrq";

    string str, odir(""), fni, strtu("s") ; //fnoc[ncmax]={"oc"}, fnob[ncmax][nbmax]={"ob"}
            // fill_n(fnoc, ncmax, "oc");
    // vector<vector<string>> fnob2 = vector<vector<string>>(ncmax,vector<string>(nbmax,"ob"));

        // string str, odir(""), fni,fnoc[ncmax]={"oc"}, fnob[ncmax][nbmax]={"ob"}, strtu("s") ;
            // fill_n(fnoc, ncmax, "oc");
            // // fill(&fnob[0][0], &fnob[ncmax][nbmax], "ob");
            // fill(&fnob[0][0], &fnob[ncmax-1][nbmax], "ob");

    double time, dt, dtminlim=1e-6, dta, dtout, outt, tend, rtu=1.0, 
        totmw=0.0, totmn=0.0, totmni=0.0, tote=0.0, 
        totmwio=0.0, totmnio=0.0, totmniio=0.0, toteio=0.0, 
        pstatnext;


    // ! imbrq : mbrq is considered if set 1, not if 0 (default)
    // ! iohsync : 1: synchronize history output of mbrq with can (default), 0: not 
    // ! ombrq : mbrq output directory (set under odir later)
    int is4, imbrq=0, iohsync=1;
    string ombrq="mbrq";
    // logical :: fex
    // ! parameters for decay heat model for input
    double t0sdown=0.0, powopr=0.0, coremass=0.0, fqdcex=0.0;
    // ! for date and time
    string dattim;

    /*##################################################################*/ 

    /*##################################################################*/ 
    /* FILE READING 67-73 */
    // namelist = READ FILE
    // & bd, cn
    #include "opr1k.dat"
    // #include "apr1k.dat"

    // ifstream input;
    
    // ifstream fin("opr1k.dat");
    /*##################################################################*/ 

    /*##################################################################*/ 
    // 76 ! initial setting for wrsttab
    wrsttabset();
    /*##################################################################*/ 


    /*##################################################################*/ 
    // 194 ! open output files
    ofstream uoc, uob;
    uoc.precision(4);
    uob.precision(4);

    
    /*##################################################################*/ 






    /*##################################################################*/ 
    // 217 ! when dt and dtout are given by time tables
    if (tdt.n > 0) {ttabcheck(tdt);}
    if (tdtout.n > 0) {ttabcheck(tdtout);}
    /*##################################################################*/ 

    /*##################################################################*/ 
    // 221 ! time unit and string
    if(itu == 1){
        strtu = "h"; 
        rtu = 1/3600;
    }
    else if(itu == 2){
        strtu = "day"; 
        rtu = 1/3600/24;
    }
    time = time/rtu; 
    tend = tend/rtu; // ! time unit conv. (input => calc.)
    /*##################################################################*/ 










    /*##################################################################*/ 
    // 229 ! write headers of output files

    // write(uoc(i),fmt='(a)') '#units:time('//trim(strtu)//'), others are in SI'
    // write(uoc(i),'(3a,i1)') '#input_file:', trim(fni), ' can #', i
    uoc.open("uoc",ios::out);
    uoc << 
            "#1:time(//trim(strtu)//) 2:pv 3:pn 4:pni 5:ptot 6:ts 7:tl 8:tlc "
            "9:volg 10:voll 11:vollc 12:rhov 13:rhon 14:rhoni 15:rhol 16:rholc "
            "17:ev 18:en 19:eni 20:el 21:elc 22:dmwev 23:mdiv 24:mdin 25:mdini 26:mdil0 "
            "27:mdil 28:mdilc 29:mdlex 30:mdov 31:mdon 32:mdoni 33:mdol 34:mdolc 35:uog "
            "36:uol 37:uolc 38:prc 39:filc 40:qsrc 41:hlc 42:hl "
            "43:totmw 44:totmn 45:totmni 46:tote 47:totmwio 48:totmnio "
            "49:totmniio 50:toteio "
            "51:totmlossv 52:totmlossn 53:totmlossni 54:totmlossl 55:totmlosslc 56:toteloss "
            "57:mdorcic" << endl;
    uoc << "\t" << fixed << time /* time*rtu */ << "\t" <<cn.pv << "\t" <<cn.pn << "\t" <<
        cn.pni << "\t" <<cn.ptot << "\t" <<cn.ts << "\t" <<cn.tl << "\t" <<cn.tlc << "\t" <<
        cn.volg << "\t" <<cn.voll << "\t" <<cn.vollc << "\t" <<cn.rhov << "\t" <<cn.rhon << "\t" <<
        cn.rhoni << "\t" << cn.rhol << "\t" <<cn.rholc << "\t" << cn.ev << "\t" <<cn.en << "\t" <<
        cn.eni << "\t" <<cn.el << "\t" <<cn.elc << "\t" <<cn.dmwev << "\t" << cn.mdiv << "\t" <<
        cn.mdin << "\t" <<cn.mdini << "\t" <<cn.mdil0 << "\t" << cn.mdil << "\t" <<cn.mdilc << "\t" <<
        cn.mdlex << "\t" <<cn.mdov << "\t" <<cn.mdon << "\t" <<cn.mdoni << "\t" << cn.mdol << "\t" <<
        cn.mdolc << "\t" <<cn.uog << "\t" << cn.uol << "\t" <<cn.uolc << "\t" <<cn.prc << "\t" <<
        cn.filc << "\t" <<cn.qsrc << "\t" <<cn.hlc << "\t" <<cn.hl << "\t" <<
        totmw << "\t" <<totmn << "\t" <<totmni << "\t" <<tote << "\t" <<totmwio << "\t" <<totmnio << "\t" <<totmniio << "\t" <<toteio << "\t" <<
        cn.totmlossv << "\t" <<cn.totmlossn << "\t" <<cn.totmlossni << "\t" <<cn.totmlossl << "\t" <<
        cn.totmlosslc << "\t" << cn.toteloss << "\t" <<cn.mdorcic << endl;
    uoc.close();


    // uob open
    // uob.open("uob",ios::out);
    // uob << 
    //         "#1:time(//trim(strtu)//) 2:fwet 3:twet 4:tdry 5:twsf 6:tdsf 7:qgen 8:qtr "
    //         "9:qrlx 10:qflxwet 11:qflxdry 12:ene 13:intqgen "
    //         "14:rho 15:cp 16:lam 17:epsrad 18:htmodw 19:htmoddw "
    //         "20:hlow 21:hhigh 22:acr 23:asf 24:tspt" << endl;
    // uob << "\t" << fixed << time /* *rtu*/<<"\t"<<bd[j-1].fwet<<"\t"<<bd[j-1].twet<<"\t"<<bd[j-1].tdry<<"\t"<<bd[j-1].twsf<<"\t"<<bd[j-1].tdsf<<"\t"<< 
    //       bd[j-1].qgen<<"\t"<<bd[j-1].qtr<<"\t"<<bd[j-1].qrlx<<"\t"<<bd[j-1].qflxwet<<"\t"<<bd[j-1].qflxdry<<"\t"<< 
    //       bd[j-1].ene<<"\t"<<bd[j-1].intqgen<<"\t"<< bd[j-1].rho<<"\t"<<bd[j-1].cp<<"\t"<<bd[j-1].lam<<"\t"<<bd[j-1].epsrad<<"\t"<< 
    //       bd[j-1].htmodw<<"\t"<< bd[j-1].htmoddw<<"\t"<< 
    //       bd[j-1].hlow<<"\t"<< bd[j-1].hhigh<<"\t"<< bd[j-1].acr<<"\t"<< bd[j-1].asf<<"\t"<< 
    //       bd[j-1].clp.tspt << endl;
    // uob.close();





    /*##################################################################*/ 








    /*##################################################################*/ 
    // 267 ! mbrq open output files and write headers
    /*##################################################################*/ 



    /*##################################################################*/ 
    // // 273 ! initialization of the can
    /* """ ONE CAN """ */
    initcan(cn);
    if(cn.vol<0.0 || cn.volg<0.0 || cn.voll<0.0 || cn.vollc<0.0){
        printf("TMP-main-initcan\n");
        // print("main: invalid volume spec. can //", 1)
        // print("  vol,volg,voll,vollc=",cn.vol,cn.volg,cn.voll,cn.vollc)
        exit(EXIT_FAILURE);
    }
    /*##################################################################*/ 










    /*##################################################################*/ 
    // 284 ! initial setting the bodies set for moving between cans
        // w: nmv = 0... 
    /*##################################################################*/ 









    /*##################################################################*/ 
    // 289 ! initialization of hot bodies
    for (i=1; i<=ncan; i++)//  do i=1, ncan
    {
        for (j=1; j<=nbody[i-1];j++) // while j < nbody: // do j=1, nbody(i) j+=1
        {
            if(bd[j-1].hlow <= 0.0 && bd[j-1].thlow.n > 0){
                ttabeval(bd[j-1].thlow, time, bd[j-1].hlow);
            }
            if(bd[j-1].hhigh <= 0.0 && bd[j-1].thhigh.n > 0){
                ttabeval(bd[j-1].thhigh, time, bd[j-1].hhigh);
            }
            if(bd[j-1].asf <= 0.0 && bd[j-1].tasf.n > 0){
                ttabeval(bd[j-1].tasf, time, bd[j-1].asf);
            }
            initbody(bd[j-1]);
        }
    }
    /*##################################################################*/ 




    /*##################################################################*/ 
    // 302 ! mbrq initialization
    if(imbrq > 0){
        ttabcheck(md.vjin); // ! check of time table inputs
        ttabcheck(md.djin);
        ttabcheck(md.tjin);
        i = md.ic;
        ifcandbrc(md,cn,1); // ! interface can in initialization mode
        initdbrc(md);

        // ! decay heat model initialization
        initdecayheat(t0sdown, powopr, coremass, fqdcex);
    }
    // end if
    /*##################################################################*/ 








    /*##################################################################*/ 
    // 315 ! debug out
    /*##################################################################*/ 










    /*##################################################################*/ 
    // 376 
    ngood = 0; // ! good step count
    ioutcount = 0; // ! output count (in case too many steps is taken for outdt)
    outt = time;
    if(tdt.n>0) {ttabeval(tdt, time, dt);} // ! if dt is given by time table, override dt
    dta = dt;
    ist = 0;
    /*##################################################################*/ 





    // TMP
    
    outdbrcool("i",ombrq, md);





















    /*##################################################################*/ 
    /*##################################################################*/ 
    // DO LOOP 382 - 614 :: if (time > tend) exit;
    // do 
    while (time < tend)
    {
        // ! if dt or dtout is given by time table, set them here
        // ! it overrides single value data dt, dtout
        if(tdt.n>0) {ttabeval(tdt, time, dt);}
        dta = min(dta, dt);
        // ! increase of dta is handled when ngood is big enough
        if(tdtout.n>0) {ttabeval(tdtout, time, dtout);}

        if(imbrq > 0 && iohsync == 1) {md.dtouth = dtout;}

        // ! solution fail => reduce dt
        // ! otherwise => count ngood, increase dt if ngood is big enough
        if(ist != 0) {
        // ! retrieve the old time data
            for (i=0;i<ncan;i++)// do i=1, ncan
            {
                cn = cnbk;
                for (j=0; j<nbody[i]; j++)//do j=1, nbody(i)
                {
                    bd[j] = bdbk[j];
                }
                // end do
            }
            // end do

            // ! retrieve the old time data for mbrq model
            if(imbrq > 0) {md = mdbk;}

            // ! reduce dt (pout overshooting is not strictly watched)
            dta = dta * 0.8;
            // write(*,'(1p,2(a,e14.6))') 'dt- ', dta, ' time ', time
            printf("'dt- ', %f, ' time ', %f\n", dta, time);
            if(dta < dtminlim) {
                // write(*,'(a)') 'main: dt too small, Aborted.'
                // write(*,'(1p,2(a,e12.4))') ' dt=', dta, ' time=',time
                printf("main: dt too small, Aborted.\n");
                printf("dt = %f,  time = %f", dta, time);
                exit(EXIT_FAILURE);
            }
            ist = 0;
            ngood = 0;
        }
        else if (ngood > 5) {
            if(dta < dt) {
                // ! increase dt
                dta = min(dta*1.2, dt);
                // write(*,'(1p,2(a,e14.6))') 'dt+ ', dta, ' time ', time
                printf("dt+  %f, time   %f\n", dta, time);
            }
            ngood = 0;
        }
        // ! back up
        // do i=1, ncan
        cnbk = cn;
        for (j=0; j<nbody[0]; j++) // do j=1, nbody(i)
        {
            bdbk[j] = bd[j];
        }
        // 433 ! back up of mbrq object
        if(imbrq > 0) {
            mdbk = md;
        }










        // 438 ! do moving of hot bodies between cans










        // 444 ! do the hot bodies
        is = 0;
        for(i=0; i<ncan; i++)// do i=1, ncan
        {
            for (j=0; j<nbody[0]; j++) // do j=1, nbody[0]
            {
                // ! pass the can's coolant conditions to hb.clp (gas 1/2/3 is hard coded as h2o,air,h2)
                bd[j].clp.ng = 3;
                bd[j].clp.gname[1]="h2o"; 
                bd[j].clp.gname[2]="air"; 
                bd[j].clp.gname[3]="h2";
                bd[j].clp.pg[1]=cn.pv; bd[j].clp.pg[2]=cn.pn; bd[j].clp.pg[3]=cn.pni;
                bd[j].clp.ptot=cn.pv + cn.pn + cn.pni;
                bd[j].clp.rhog=cn.rhov + cn.rhon + cn.rhoni;
                bd[j].clp.tl=cn.tl; bd[j].clp.tg=cn.ts;
                bd[j].clp.rhol=cn.rhol; bd[j].clp.cpl=cn.cpl; bd[j].clp.betal=cn.betal;
                // ! gas phys. props. need consideration on mixture, not given here
                bd[j].clp.hlev=cn.hl;
                bd[j].clp.wvol=cn.voll; bd[j].clp.wmdi=cn.mdil;
                // ! do the evolution
               hbevol(bd[j], time, dta, is);
                if(is != 0) {
                    ist = is; 
                    break;
                }
                // end if
            }
            if(ist != 0) {break;}
        }
        // printf("main-444\n");
        printf("time < tend : %f < %f\n", time, tend);
        if(ist != 0) {continue;} // ! if solution fails, quit here and redo with reduced dt

        // 469 ! time advancement of mbrq
        if(imbrq > 0) {
            decayheat(time); // ! evaluate decay heat
            i = md.ic;
            ifcandbrc(md,cn);   // ! interface hotcan data -> mbrq
            dbrcevol(md,time,dta); // ! evolution of mbrq
        }

        // 477 ! do the hot cans
        // ! do the hot cans
        // ! ONE CAN
        // for (i=1; i<=ncan; i++) // do i=1, ncan
        // {
        // ! set the heat source for the can
        cn.qsrcex = 0.0;
        cn.tsrcex = cn.ts;
        for (j=1; j<=nbody[0]; j++) // do j=1, nbody[i-1]
        {
            cn.qsrcex = cn.qsrcex + bd[j-1].qtr;
            // ! setting the heat source temperature for checking temp. overshooting
            if(cn.qsrcex > 0.0) {
                if(bd[j-1].fwet < 1.0) {
                    cn.tsrcex = max(cn.tsrcex, bd[j-1].tdsf);
                }
                else {
                    cn.tsrcex = max(cn.tsrcex, bd[j-1].twsf);
                }
            }
            else {
                if(bd[j-1].fwet < 1.0) {
                    cn.tsrcex = min(cn.tsrcex, bd[j-1].tdsf);
                }
                else {
                    cn.tsrcex = min(cn.tsrcex, bd[j-1].twsf);
                }
            }
        }

        // 500 ! heat input from mbrq
        if(imbrq > 0){
            if(1 == md.ic ) {ifdbrccan(md,cn,time, dta);} // ! interface mbrq -> hotcan
            md.totqsrcex = md.totqsrcex + cn.qsrcex * dta; // ! counting total external heat source
        } // ! (mbrq + body)

        // ONE CAN // 
        // ! set the can to can transfers
        // ! the exit mass flow (calculated) is fed to the downstream can directly.
        // ! the inlet flow given by the input is alive, the transfer from the
        // ! upstream can is added to them.
        // ! so, the inlet mass flow ourput (mdi*) does not show what is given here.
        // if(i < ncan){ // ! connection of the pressures
        //     preslevel(cn(i+1),cn(i).poutlev, pstatnext) 
        //     cn(i).pout = pstatnext
        // end if
        // if(i > 1){ // ! receive flow-out, work-out from the up-stream can
        //     cn(i).mv = cn(i).mv + cn(i-1).mdov*dta
        //     cn(i).mn = cn(i).mn + cn(i-1).mdon*dta
        //     cn(i).mni = cn(i).mni + cn(i-1).mdoni*dta
        //     cn(i).ml = cn(i).ml + cn(i-1).mdol*dta
        //     cn(i).mlc = cn(i).mlc + cn(i-1).mdolc*dta
        //     cn(i).etot = cn(i).etot + &
        //                 (cn(i-1).mdov*cn(i-1).ev + cn(i-1).mdon*cn(i-1).en &
        //             + cn(i-1).mdoni*cn(i-1).eni &
        //             + cn(i-1).mdol*cn(i-1).el + cn(i-1).mdolc*cn(i-1).elc &
        //             + cn(i-1).wko ) * dta
        // end if       
        // if(i < ncan .and. cn(i+1).tmdorcic.n > 0){ // ! RCIC connection (revers flow)
        //     cn(i).ml = cn(i).ml + cn(i+1).mdorcic * dta
        //     cn(i).etot = cn(i).etot + cn(i+1).mdorcic * cn(i+1).el* dta
        // end if
        // end do // ! untill here, all the process is about the old time values

        // 533 ! do the update of cans now
        // do i=1, ncan

        massenevol(cn, time, dta, xio);

        // if (time > 65.8){
        //     printf("here\n");
        // }


        getequil(cn, 0, is);
        // printf ("main-here\n");

        if(is != 0 && !(is == -2 && dta<1.0)) { // ! pout overshooting with
            ist = is; 
            break;       // ! successful solution allowed
        }
        // end do
        if(ist != 0) {continue;} // ! if solution fails, quit here and redo with reduced dt




        // 543 ! do the following when the solution is success
        ngood = ngood + 1; // ! count the good cycle

        for (i=1; i<=ncan; i++) // do i = 1, ncan
        {
            totmw = cn.mv + cn.ml + cn.mlc;
            totmn = cn.mn ; 
            totmni = cn.mni;
            tote = cn.mv*cn.ev + cn.ml*cn.el + cn.mlc*cn.elc 
                        + cn.mn*cn.en + cn.mni*cn.eni;
            totmwio = totmwio + xio.mwio;
            totmnio = totmnio + xio.mnio;
            totmniio = totmniio + xio.mniio;
            toteio = toteio + xio.eio;
        }

        time = time + dta;
        istep = istep + 1;
        ioutcount = ioutcount + 1;











































        if( time >= outt || ioutcount > ostep || time > tend){
            printf("WRITE TO THE FILE:: oc & ob \n");
            // if(icmma == 1){
            //     str='(f15.6,1p,56(",",e14.6e3))'
            // else
            //     str='(f15.6,1p,56(1x,e14.6e3))'
            // end if
            // do i = 1, ncan
                uoc.open("uoc",ios::app);
                // write(uoc[i-1],fmt=str) 
                uoc << "\t" << fixed << time /* time*rtu */ << "\t" <<cn.pv << "\t" <<cn.pn << "\t" <<
                    cn.pni << "\t" <<cn.ptot << "\t" <<cn.ts << "\t" <<cn.tl << "\t" <<cn.tlc << "\t" <<
                    cn.volg << "\t" <<cn.voll << "\t" <<cn.vollc << "\t" <<cn.rhov << "\t" <<cn.rhon << "\t" <<
                    cn.rhoni << "\t" << cn.rhol << "\t" <<cn.rholc << "\t" << cn.ev << "\t" <<cn.en << "\t" <<
                    cn.eni << "\t" <<cn.el << "\t" <<cn.elc << "\t" <<cn.dmwev << "\t" << cn.mdiv << "\t" <<
                    cn.mdin << "\t" <<cn.mdini << "\t" <<cn.mdil0 << "\t" << cn.mdil << "\t" <<cn.mdilc << "\t" <<
                    cn.mdlex << "\t" <<cn.mdov << "\t" <<cn.mdon << "\t" <<cn.mdoni << "\t" << cn.mdol << "\t" <<
                    cn.mdolc << "\t" <<cn.uog << "\t" << cn.uol << "\t" <<cn.uolc << "\t" <<cn.prc << "\t" <<
                    cn.filc << "\t" <<cn.qsrc << "\t" <<cn.hlc << "\t" <<cn.hl << "\t" <<
                    totmw << "\t" <<totmn << "\t" <<totmni << "\t" <<tote << "\t" <<totmwio << "\t" <<totmnio << "\t" <<totmniio << "\t" <<toteio << "\t" <<                    
                    cn.totmlossv << "\t" <<cn.totmlossn << "\t" <<cn.totmlossni << "\t" <<cn.totmlossl << "\t" <<
                    cn.totmlosslc << "\t" << cn.toteloss << "\t" <<cn.mdorcic << endl;
                uoc.close();
            // end do

            // if(icmma == 1){
            //     str='(f15.6,1p,16(",",e14.6e3),2(",",a),5(",",e14.6e3))'
            // else
            //     str='(f15.6,1p,16(1x,e14.6e3),2(1x,a),5(1x,e14.6e3))'
            // end if
            // do i=1, ncan
            // do j=1, nbody[i-1]
            for (j=1; j<=nbody[0]; j++){

            //     write(uob(i,j),fmt=str) 
            //         time*rtu,bd[j-1].fwet,bd[j-1].twet,bd[j-1].tdry,bd[j-1].twsf,bd[j-1].tdsf, 
            //         bd[j-1].qgen,bd[j-1].qtr,bd[j-1].qrlx,bd[j-1].qflxwet,bd[j-1].qflxdry, 
            //         bd[j-1].ene,bd[j-1].intqgen, bd[j-1].rho,bd[j-1].cp,bd[j-1].lam,bd[j-1].epsrad, 
            //         bd[j-1].htmodw, bd[j-1].htmoddw, 
            //         bd[j-1].hlow, bd[j-1].hhigh, bd[j-1].acr, bd[j-1].asf, 
            //         bd[j-1].clp.tspt 
            }
            // end do
            // end do

            // outt = outt + dtout
            // ioutcount = 0
        }

        // 603 ! mbrq output
        if(imbrq > 0) {
            if(time > tend) { // ! force write history data at the end step
                outdbrcool("W",ombrq,md, time, istep);
            }
            else { // ! otherwise the proper write timing is checked
                outdbrcool("w",ombrq,md, time, istep);
            }
        }


    }   // 612
    /*##################################################################*/ 
    /*##################################################################*/ 

    // 616 ! close files

    if (imbrq > 0) {
        outdbrcool("c",ombrq, md);
        outparsieve(ombrq,md);
    }






    /*##################################################################*/ 

    printf("hello1 \n");
    printf("hello2 \n");

    // // melthex mh;
    // // mh.nthex = 3;
    // // float time=1, dt=0.1;
    // // std::cout << typeid(dt).name() << std::endl;
    // md.psieve.nn=20;
    // printf ("\n %d \n", md.psieve.nn);
    // // mjetevol(&md, time, dt);
    // mjetevol(md, time, dt);
    // printf ("\n %d \n", md.psieve.nn);
    

    // printf(" ################\n ");
    // printf("DONE\n");

    return 0;
}