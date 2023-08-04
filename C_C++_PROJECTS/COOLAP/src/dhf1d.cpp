#ifndef DHF1D_H
#define DHF1D_H

#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

#include "consts.cpp"
#include "wrsttab.cpp"
// #include "misc.cpp"

/*###############################################################################*/ 
// ! calculate 1D debris bed DHF
// ! need to initialize steam table by "wrsttab" at the beginning
/*###############################################################################*/ 
double pressure, temp, porosity, dp, 
    sigma, den_l, den_g, mu_l, mu_g, K_per, n_pas, del_h, wr, 
    pcr=22.06e6;
struct physprop_dhf
{
    int ip = 1;
        // ! 0: water at subcool/superheat
        // ! 1: water at saturation at given pressure (default)
        // ! 2: water at saturation at given temp
        // ! 3: given by user
    double rhol,rhog,mul,mug,dhlg,sig;
        // ! these are used only for ip=3, 
        // ! otherwise water prop. for p,t are calculated
};
/*###############################################################################*/ 

/*###############################################################################*/ 
void setprop_dhf(double &p,double &t,double por,double diap,double vl)//, phpr)
{
    double rl,rg,drdp,drdt,hl,hg,dhdp,dhdt,s;
    int irng, ip = 1; // ! phpr%ip default=sat. given p;
    
    //   if(present(phpr)) ip = phpr%ip
    if(ip == 0) { // ! subcool/superheat
        wrsteamtab(p,t,0, rl,drdp,drdt, hl,dhdp,dhdt,irng);
        wrsteamtab(p,t,1, rg,drdp,drdt, hg,dhdp,dhdt,irng);
        den_l = rl; 
        den_g = rg; 
        del_h = hg - hl;
        mu_l = viscwaterf(rl,t); 
        mu_g = viscwaterf(rg,t);
        sigma = stenwaterf(t);
    }
    else if(ip == 1){ // ! sat. given p
        t = tsatwaterf(p);
        wrsteamtab(p,t,0, rl,drdp,drdt, hl,dhdp,dhdt,irng);
        wrsteamtab(p,t,1, rg,drdp,drdt, hg,dhdp,dhdt,irng);
        den_l = rl; 
        den_g = rg; 
        del_h = hg - hl;
        mu_l = viscwaterf(rl,t); 
        mu_g = viscwaterf(rg,t);
        sigma = stenwaterf(t);
    }
    else if(ip == 2){ // ! sat. given T
        p = psatwaterf(t);
        wrsteamtab(p,t,0, rl,drdp,drdt, hl,dhdp,dhdt,irng);
        wrsteamtab(p,t,1, rg,drdp,drdt, hg,dhdp,dhdt,irng);
        den_l = rl; 
        den_g = rg; 
        del_h = hg - hl;
        mu_l = viscwaterf(rl,t); 
        mu_g = viscwaterf(rg,t);
        sigma = stenwaterf(t);
    }
    else if(ip == 3){ // ! arbitrarily given properties
        printf(" # No phpr - whang \n");
        exit(EXIT_FAILURE);
        // sigma = phpr.sig;
        // den_l = phpr.rhol; den_g = phpr.rhog;
        // mu_l = phpr.mul; mu_g = phpr.mug;
        // del_h = phpr.dhlg;
    }
    else {
        printf(" # wrong phys.prop option\n");
        exit(EXIT_FAILURE);
        // write(*,'(a)') '# wrong phys.prop option'
        // stop
    }
    pressure = p;
    temp = t;
    porosity = por ; 
    dp = diap;
    wr = vl;
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void part_gas_drag(double jg, double alpha, double &ff, int fmod)
{
    double K_g, n_g, a2, a3, a4, EA, zeta, W, dp_h, a5, a6;
    int fn = 3; // ! default model: Rahman

    fn = fmod;
    // ! REED, Lip, MTD_Rahman, MTD_Lee
    if (fn==1){ // ! REED
        K_g=pow(alpha,3);
        n_g=pow(alpha,5);
        ff=alpha*porosity*(mu_g/K_per/K_g*jg+den_g/n_pas/n_g*fabs(jg)*jg);
    }
    else if (fn==2){ // ! Lip
        K_g=pow(alpha,3);
        n_g=pow(alpha,3);
        ff=alpha*porosity*(mu_g/K_per/K_g*jg+den_g/n_pas/n_g*fabs(jg)*jg);
    }
    else if (fn==3){ // ! MTD_Rahman
        a3=max(0.0,min(0.6,0.6+(pow((dp-0.012),3))*4e5));
        a4=max(0.0,min(0.74,(pow((dp-0.012),3))*4e5+0.74));
        EA=(1.0-porosity)/(1.0-porosity*alpha);
        zeta=(alpha-a3)/(a4-a3);
        W=pow(zeta,2)*(3.0-2.0*zeta);
        if (alpha<=a3){ 
            K_g=pow(alpha,4);
            n_g=pow(alpha,4);
        }
        else if (alpha<=a4){ 
            K_g=pow(alpha,4)/(1.0-W*(1.0-alpha));
            n_g=pow(alpha,4)/(1.0-W*(1.0-alpha));
        }
        else {
            K_g=pow(alpha,3);
            n_g=pow(alpha,3);
        }
        ff=porosity*alpha*(mu_g/K_per/K_g*jg+den_g/n_pas/n_g*fabs(jg)*jg);
    }
    else if (fn==4){ // ! MTD_Lee
        dp_h=1.5*porosity/(1.0-porosity)*dp;
        a3=max(0.0,min(0.6,0.6+(pow((dp_h-0.012),3))*4e5));
        a4=max(0.0,min(0.74,(pow((dp_h-0.012),3))*4e5+0.74));
        a5=max(0.0,min(0.8,(dp_h*1e3-7.0)/2.0e1+0.8));
        a6=1.0;
        EA=(1.0-porosity)/(1.0-porosity*alpha);

        if (alpha<=a3){ 
            K_g=pow(alpha,4);
            n_g=pow(alpha,4);
        }
        else if (alpha<=a4){ 
            zeta=(alpha-a3)/(a4-a3);
            W=pow(zeta,2)*(3.0-2.0*zeta);
            K_g=pow(alpha,4)/(1.0-W*(1.0-alpha));
            n_g=pow(alpha,4)/(1.0-W*(1.0-alpha));
        }
        else if (alpha<=a5){ 
            K_g=pow(alpha,3);
            n_g=pow(alpha,3);
        }
        else {
            zeta=(alpha-a5)/(a6-a5);
            W=pow(zeta,2)*(3.0-2.0*zeta);
            K_g=((pow(EA,(4./3.)))*pow(alpha,3))/((pow(EA,(4./3.)))+W*(1.0-(pow(EA,(4./3.)))));
            n_g=((pow(EA,(2./3.)))*pow(alpha,3))/((pow(EA,(2./3.)))+W*(1.0-(pow(EA,(2./3.)))));
        }
        ff=porosity*alpha*(mu_g/K_per/K_g*jg+den_g/n_pas/n_g*fabs(jg)*jg);
    }
    else {
        // write(*,*) 'part_gas_drag: unkonwn model: fn=',fn
        // stop
        printf("part_gas_drag: unkonwn model: fn=%d\n",fn);
        exit(EXIT_FAILURE);
    }
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void part_liq_drag(double jl, double alpha, double &ff, int fmod)
{
    double K_l, n_l;
    int fn = 3; // ! default model: Rahman

    fn = fmod;
    // ! REED, Lip, MTD_Rahman
    if (fn==1){ // ! REED
        K_l=pow((1.0-alpha),3);
        n_l=pow((1.0-alpha),5);
        ff=(1.0-alpha)*porosity*(mu_l/K_per/K_l*jl+den_l/n_pas/n_l*fabs(jl)*jl);
    }
    else if (fn==2){ // ! Lip
        K_l=pow((1.0-alpha),3);
        n_l=pow((1.0-alpha),3);
        ff=(1.0-alpha)*porosity*(mu_l/K_per/K_l*jl+den_l/n_pas/n_l*fabs(jl)*jl);
    }
    else if (fn==3){ // ! MTD_Rahman
        K_l=pow((1.0-alpha),3);
        n_l=pow((1.0-alpha),6);
        ff=porosity*(1.0-alpha)*(mu_l/K_per/K_l*jl+den_l/n_pas/n_l*fabs(jl)*jl);
    }
    else if (fn==4){ // ! MTD_Lee
        K_l=pow((1.0-alpha),4);
        n_l=pow((1.0-alpha),4);
        ff=(1.0-alpha)*porosity*(mu_l/K_per/K_l*jl+den_l/n_pas/n_l*fabs(jl)*jl);
    }
    else {
        // write(*,*) 'part_liq_drag: unkonwn model: fn=',fn
        // stop
        printf ("part_liq_drag: unkonwn model: fn= %d \n", fn);
        exit(EXIT_FAILURE);
    }
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void glint_drag(double jl,double jg, double alpha, double &ff, int fmod)
{
    double coeff_force, eff_vel_g, eff_vel_l, s, den_m, js, f,
        a0,a1,a2,a3,a4, ratio, Db, eta, jr, EA, a_bar, b_bar, 
        Ca_annular, Cb_annular, C_v, C_i, Ca, Cb, zeta2, W2, zeta3, W3, mp,
        dp_h, a_bar_ann, b_bar_ann, Ca_channel, Cb_channel, zeta4, W4, 
        FNA,FNB, a5, a6;
    int fn = 3; // ! default model: Rahman

    fn = fmod;

    // ! REED, Lip, MTD_Rahman
    if (fn==1 || fn==2){ // ! REED/Lip no g-l interaction
        ff = 0.0;
    }
    else if (fn==3){  // ! MTD_Rahman
        Db = 1.35*(pow((sigma/fabs(grav)/(den_l-den_g)),0.5));
        ratio = Db/dp;
        eta = pow((pi*(pow(2.0,0.5))/6.0/(1.0-porosity)),(1./3.));

        a0=max((pi/3.0)*((1.0-porosity)/porosity)*ratio*(ratio+1.0)* 
            (6.0*eta-5.0*(1.0+ratio)), 0.0);
        a1=max(0.0,min(0.3,0.3+(pow((dp-0.012),3))*4e5));
        a2=max(0.0,min(pi/6.0,pi/6.0+(pow((dp-0.012),3))*4e5));
        a3=max(0.0,min(0.6,0.6+(pow((dp-0.012),3))*4e5));
        a4=max(0.0,min(0.74,(pow((dp-0.012),3))*4e5+0.74));

        den_m = den_l*((1.0-alpha)+alpha*den_g/den_l);
        js=jg/alpha-jg-jl;
        f=0.5*(1.0+ratio)*log(1.0+2.0/ratio);

        Ca= (mu_l/Db/Db)*min(1.0,dp/0.012);
        Cb= ((den_l*(1.0-alpha)+den_g*alpha)/Db/porosity)*min(1.0,dp/0.012);

        mp=0.25*min(1.0,pow((dp/0.003),3));

        EA=(1.0-porosity)/(1.0-porosity*alpha);

        jr=jg-alpha/(1.0-alpha)*jl;
        a_bar = mp/(K_per*(pow(alpha,3))*(pow(1.0,(4./3.))));
        b_bar = mp/(n_pas*(pow(alpha,3))*(pow(1.0,(2./3.))));

        if (alpha<= a0){ 
            C_v=18.0*alpha*f;
            C_i=0.34*pow((1.0-alpha),3)*alpha*pow(f,2);
            ff=Ca*C_v*js+Cb*C_i*fabs(js)*js;
        }
        else if (alpha<= a1){ 
            C_v=18.0*(a0*f+alpha-a0);
            C_i=0.34*pow((1.0-alpha),3)*(a0*pow(f,2)+alpha-a0);
            ff=Ca*C_v*js+Cb*C_i*fabs(js)*js;
        }
        else if (alpha<= a2){ 
            zeta2=(alpha-a1)/(a2-a1);
            W2=pow(zeta2,2)*(3.0-2.0*zeta2);
            C_v=18.0*(a0*f+alpha-a0)*(1.0-W2)+5.21*alpha*W2;
            C_i=0.34*pow((1.0-alpha),3)*(a0*pow(f,2)+alpha-a0)*(1.0-W2)+ 
                W2*0.92*pow((1.0-alpha),3)*alpha;
            ff=Ca*C_v*js+Cb*C_i*fabs(js)*js;
        }
        else if (alpha<= a3){ 
            C_v=5.21*alpha;
            C_i=0.92*pow((1.0-alpha),3)*alpha;
            ff=Ca*C_v*js+Cb*C_i*fabs(js)*js;
        }
        else if (alpha<= a4){ 
            zeta3=(alpha-a3)/(a4-a3);
            W3=pow(zeta3,2)*(3.0-2.0*zeta3);
            C_v=5.21*alpha;
            C_i=0.92*pow((1.0-alpha),3)*alpha;
            Ca_annular= porosity*(1.0-alpha)*mu_g*a_bar;
            Cb_annular= porosity*(1.0-alpha)*den_g*b_bar;
            ff=Ca*C_v*(1.0-W3)*js+W3*Ca_annular*jr+ 
                (Cb*C_i*(1.0-W3)*fabs(js)*js+W3*Cb_annular*fabs(jr)*jr);
        }
        else if (alpha<=1.0){ 
            ff=porosity*(1.0-alpha)*(mu_g*a_bar+fabs(jr)*den_g*b_bar)*jr;
        }
        else {
            ff=0.0;
        }
    }
    else if (fn==4){  // ! MTD_Lee
        dp_h=1.5*porosity/(1.0-porosity)*dp;
        Db = min(1.35*(pow((sigma/fabs(grav)/(den_l-den_g)),0.5)),0.41*dp_h);
        ratio = Db/dp_h;
        eta = pow((pi*(pow(2.0,0.5))/6.0/(1.0-porosity)),(1./3.));

        a0=max((pi/3.0)*((1.0-porosity)/porosity)*ratio*(ratio+1.0)* 
            (6.0*eta-5.0*(1.0+ratio)), 0.0);
        a1=max(0.0,min(0.3,0.3+(pow((dp_h-0.012),3))*4e5));
        a2=max(0.0,min(pi/6.0,pi/6.0+(pow((dp_h-0.012),3))*4e5));
        a3=max(0.0,min(0.6,0.6+(pow((dp_h-0.012),3))*4e5));
        a4=max(0.0,min(0.74,(pow((dp_h-0.012),3))*4e5+0.74));
        a5=max(0.0,min(0.8,(dp_h*1e3-7.0)/2.0e1+0.8));
        a6=1.0;

        den_m = den_l*((1.0-alpha)+alpha*den_g/den_l);
        js=jg/alpha-jg-jl;
        f=0.5*(1.0+ratio)*log(1.0+2.0/ratio);

        Ca= (mu_l/Db/Db)*min(1.0,dp_h/0.012);
        Cb= ((den_l*(1.0-alpha)+den_g*alpha)/Db/porosity)*min(1.0,dp_h/0.012);

        EA=(1.0-porosity)/(1.0-porosity*alpha);

        jr=jg-alpha/(1.0-alpha)*jl;

        a_bar = (1.0/(K_per*(pow(alpha,3))))*(pow((1.0-alpha),2))*min(1.0,(dp_h/0.008));
        b_bar = (1.0/(n_pas*(pow(alpha,3))))*(pow((1.0-alpha),2))*min(1.0,(dp_h/0.008));

        a_bar_ann = (1.0/(K_per*(pow(alpha,3))*(pow(EA,(4./3.)))))*min(1.0,(dp_h/0.008));
        b_bar_ann = (1.0/(n_pas*(pow(alpha,3))*(pow(EA,(2./3.)))))*min(1.0,(dp_h/0.008));

        if (alpha<= a0){ 
            C_v=18.0*alpha*f;
            C_i=0.34*pow((1.0-alpha),3)*alpha*pow(f,2);
            ff=Ca*C_v*js+Cb*C_i*fabs(js)*js;
        }
        else if (alpha<= a1){ 
            C_v=18.0*(a0*f+alpha-a0);
            C_i=0.34*pow((1.0-alpha),3)*(a0*pow(f,2)+alpha-a0);
            ff=Ca*C_v*js+Cb*C_i*fabs(js)*js;
        }
        else if (alpha<= a2){ 
            zeta2=(alpha-a1)/(a2-a1);
            W2=pow(zeta2,2)*(3.0-2.0*zeta2);
            C_v=18.0*(a0*f+alpha-a0)*(1.0-W2)+5.21*alpha*W2;
            C_i=0.34*pow((1.0-alpha),3)*(a0*pow(f,2)+alpha-a0)*(1.0-W2)+ 
                W2*0.92*pow((1.0-alpha),3)*alpha;
            ff=Ca*C_v*js+Cb*C_i*fabs(js)*js;
        }
        else if (alpha<= a3){ 
            C_v=5.21*alpha;
            C_i=0.92*pow((1.0-alpha),3)*alpha;
            ff=Ca*C_v*js+Cb*C_i*fabs(js)*js;
        }
        else if (alpha<= a4){ 
            zeta3=(alpha-a3)/(a4-a3);
            W3=pow(zeta3,2)*(3.0-2.0*zeta3);
            C_v=5.21*alpha;
            C_i=0.92*pow((1.0-alpha),3)*alpha;
            Ca_channel= porosity*(1.0-alpha)*mu_g*a_bar;
            Cb_channel= porosity*(1.0-alpha)*den_g*b_bar;
            ff=Ca*C_v*(1.0-W3)*js+W3*Ca_channel*jr+ 
                (Cb*C_i*(1.0-W3)*fabs(js)*js+W3*Cb_channel*fabs(jr)*jr);
        }
        else if (alpha<=a5){ 
            ff=porosity*(1.0-alpha)*(mu_g*a_bar+fabs(jr)*den_g*b_bar)*jr;
        }
        else if (alpha<= a6){ 
            zeta4=(alpha-a5)/(a6-a5);
            W4=pow(zeta4,2)*(3.0-2.0*zeta4);
            FNA=porosity*(1.0-alpha)*(mu_g*a_bar+fabs(jr)*den_g*b_bar)*jr;
            FNB=porosity*(1.0-alpha)*(mu_g*a_bar_ann+fabs(jr)*den_g*b_bar_ann)*jr;
            ff=FNA*(1.0-W4)+FNB*W4;
            }
        else {
            ff=0.0;
        }
    }
    else {
        // write(*,*) 'glint_drag: unkonwn model: fn=',fn
        // stop
        printf("glint_drag: unkonwn model: fn= %d\n",fn);
        exit(EXIT_FAILURE);
    }
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void balanceeq(double alp,double qflx,double &ff,double &ffres,int fmod)
{
    // ! fmod : model no.
    // ! alp : void fraction
    // ! qflx : DHF
    // ! ff : value of balance eq
    // ! dff : dff/dalp for Newton iteration
    double x1, jj_g,jj_l, Fpga,Fpla,Fia,nrm;
    int fn = 3; // ! default model: Rahman

    fn = fmod;

    x1=fabs(grav)*(den_l-den_g)*porosity;
    jj_g = qflx/del_h/den_g;
    jj_l = wr-qflx/del_h/den_l;
    part_gas_drag(jj_g,alp, Fpga, fn);
    Fpga = Fpga/x1;
    part_liq_drag(jj_l,alp, Fpla, fn);
    Fpla = Fpla/x1;
    glint_drag(jj_l,jj_g,alp, Fia, fn);
    Fia = Fia/x1;
    ff = -alp*(1.0-alp)+(1.0-alp)*Fpga-alp*Fpla+Fia;
    nrm = fabs(alp*(1.0-alp))+fabs((1.0-alp)*Fpga)+fabs(alp*Fpla)+fabs(Fia);
    ffres = fabs(ff)/max(nrm,1e-12); // ! measure of residual of ff
    // !DEBUG
    // !  if(abs(ff) > 1d16) then
    // !    errorhdr(mess='balanceeq: ff>1d16',act='cont')
    // !  endif
}
/*###############################################################################*/ 

/*###############################################################################*/ 
void get_alp_dhf(double &alpmax, double &dhf, int &imod, int fmod, double alpini)
{
    // imod ! 1:Reed, 2:Lipinski, 3:Rahman, 4: MTD_Lee
    // alpini
    // ! initial guess alp: used if non-zero value is given
    // ! bi-section is sensitive to initial guess of alp because
    // ! existance of solution for DHF depends on it
    // ! by default, built in correlation is used

    double x1, alp, dalp, 
        qflx, qf0, qf1, ff, dff, ff0, ff1, ffres, alpa, // al(0:2), qf(0:2), 
        dqf=100e3, qfmax=1e8, alpguess;
        // ! dqf : span for scanning DHF range (W/m2)
        // ! qfmax : max reasonable limit for DHF (W/m2)
        // ! alpguess : initial guess of void frac. at DHF; solution success/fail sensitive to this
    int looplimb = 100, looplima = 100, ib, ibfail, ia, ifallback, fn, iuserguess=0;
        // ! looplimb, looplimn: limit for bi-section and Newton loop
        // ! ib, in : loop counter
        // ! ibfail : bi-section fail flag
        // ! fn : model index (default 3=Rahman, 1=Reed)
        // ! ifallback : flag of fall back to simpler Reed model (fn=1)
        // ! iuserguess ; flag for user input of alpini
    vector<double> al = vector<double>(3,0.0);
    vector<double> qf = vector<double>(3,0.0);

    fn = 3; // ! Rahman model as default
    ifallback = 0; // ! clear fallback flag

    // if(present(fmod)) fn = fmod // ! default overrided if given
    fn = fmod; 
    // if(present(alpini)){
        if(alpini > 0e0) {iuserguess = 1;} // ! a valid user alpini given
    // end if

    if(iuserguess > 0){alpguess = alpini;}
    else { // ! alp guess based on pressure dependence of DHF models
        if(fn == 3){ // ! for Rahman model
            alpguess = -0.105*log10(pressure/pcr)+0.495;
        }
        else if(fn == 1){ // ! for Reed model
            alpguess = -0.0900*log10(pressure/pcr)+0.572;
        }
        else if(fn == 4){ // ! for MTD_Lee model
            alpguess = -0.0858*log10(pressure/pcr)+0.565;
        }
        else {
            // errorhdr(mess='get_alp_dhf: unknown model',act='stop')
            printf("get_alp_dhf: unknown model\n");
            exit(EXIT_FAILURE);
        }
    }

    K_per = (pow(porosity,3))*(pow(dp,2))/(150.0*(pow((1.0-porosity),2))); // ! permeability
    n_pas = (pow(porosity,3))*dp/(1.75*(1.0-porosity)); // ! passability

    // ! Numerical convergence: Bisection for qflx, alp succesive correction

    alp = alpguess; 
    dalp = 0.2; // ! initial guess
    ia = -1;
    do                    // ! loop for alp convergence
    {
        ia = ia + 1; // ! starts with 0
        if( ia%3 == 0){  // ! create a set of 3 alp values
            al[0] = max(alp-dalp,0.001);
            al[1] = alp;
            al[2] = min(alp+dalp,0.999);
            if(al[0] >= alp) {al[0] = alp/2.0;}
            if(al[2] <= alp) {al[2] = (1.0+alp)/2.0;}
        }
        alp = al[ia%3];

        // ! setting initial two points
        qf0 = 1.0; 
        qf1 = dqf;
        do
        {
            balanceeq(alp,qf0,ff0,ffres, fn);
            balanceeq(alp,qf1,ff1,ffres, fn);

            if( (ff0*ff1) > 0.0){
                qf0 = qf1; 
                qf1=qf1+dqf;
                if(qf1 > qfmax){ // ! too large DHF, no solution available
                    if(ifallback == 1) {
                        // errorhdr(mess='get_alp_dhf: bi-sec. fail',act='stop')
                        printf("get_alp_dhf: bi-sec. fail\n");
                        exit(EXIT_FAILURE);
                    }
                    // ! fallback to Reed model
                    // errorhdr(mess='get_alp_dhf: bi-sec. try Reed model',act='cont')
                    printf("get_alp_dhf: bi-sec. try Reed model\n");
                    printf("ERRORHDR: act = cont?? \n ");
                    fn = 1; 
                    ifallback=1; 
                    qf0 = 1.0; 
                    qf1 = dqf; 
                    dalp=0.1;
                    alpguess = -0.09*log10(pressure/pcr)+0.572;
                    alp = alpguess;
                } 
            }
            else {break;}
        }
        while(ia>=0);

        ib = 0;
        ibfail = 0;
        do                 // ! loop for bisection to solve qflx
        {
            ib = ib + 1;
            qflx = (qf0 + qf1)/2.0;
            balanceeq(alp,qflx,ff,ffres, fn);

            if(ff*ff0 <= 0.0){
                qf1 = qflx; 
                ff1 = ff;
            }
            else if(ff*ff1 <= 0.0){
                qf0 = qflx; 
                ff0 = ff;
            }
            else {
                // errorhdr(mess='get_alp_dhf: bi-sec. fail (in loop)',act='stop')
                printf("get_alp_dhf: bi-sec. fail (in loop)\n");
                exit(EXIT_FAILURE);
            }
            if(ffres < 1e-7) {break;}
            if(ib >= looplimb){
                // write(*,'(a,1p,3(1x,e12.4))') 
                // 'get_alp_dhf: bisec bad convergence alp,qflx,ff=',alp,qflx,ff
                printf("get_alp_dhf: bisec bad convergence alp,qflx,ff=%f, %f, %f\n", alp,qflx,ff);
                ibfail = 1;
                break;
            }
        }
        while(ib>=0);
        // !DEBUG
        // !    write(*,'(a,i3,1x,f8.6,1x,1pe13.6)') 
        // !      'bisec: ib,alp,qflx', ib, alp, qflx
        // !END DEBUG

        qf[ia%3] = qflx;

        // ! every set of 3 alp values, seek the peak position by
        // ! quatratic function approximation of (alp,qflx)
        if(ia%3 == 2){
            // ! peak position of quadratic
            x1 = al[0]*(qf[1]-qf[2])+al[1]*(qf[2]-qf[0])+al[2]*(qf[0]-qf[1]);
            if (x1 != 0.0){
                alpa = 0.5*(pow(al[0],2)*(qf[1]-qf[2])+pow(al[1],2)*(qf[2]-qf[0])+ 
                    pow(al[2],2)*(qf[0]-qf[1])) / x1; // ! guess of alp at qf peak
            }
            else if(qf[0] == qf[1] && qf[0] == qf[2]){
                alpa = al[1]; // ! converged qf regardless of alp
            }
            else {alpa = (al[1] + 0.75)/2.0;} // ! give a new guess
            // !DEBUG
            // !      write(*,'(a,4(1x,f7.4),a,1p,2(1x,e11.4))') 
            // !        'al0,al1,al2,alpa/dqf0,dqf2',al[0],al[1],al[2],alpa,' /',
            // !        qf[0]-qf[1],qf[2]-qf[1]
            // !END DEBUG
            // ! if the peak is in the range al[0]--al[2], modify dalp according to the range
            if((qf[1]>=qf[0] && qf[1]>=qf[2]) && (alpa>=al[0] && alpa<=al[2])) 
            {
                dalp = max(fabs(alpa - al[1]),1e-5);
            }
            // ! convergence check
            if((fabs(alpa-al[1]) < 1e-4 || fabs(2.0*qf[1]-qf[0]-qf[2])/max(qf[1],1.) < 1e-4) 
                && ibfail == 0) {break;}

            // ! update alp
            alp = max(1e-4, min(0.9999, alpa));
            // ! check too slow or failed conversion
            if(ia >= looplima){
                // write(*,'(a,1p,5(1x,e12.4))') 
                // 'get_alp_dhf: alp bad convergence alp,dalp,qf(0:2)=',alp,dalp,qf[0],qf[1],qf[2]
                // exit
                printf("get_alp_dhf: alp bad convergence alp,dalp,qf(0:2)=%f %f %f %f %f\n",alp,dalp,qf[0],qf[1],qf[2]);
                break;
            }
        }
    }
    while (ia>=0);

    dhf = qflx;
    alpmax = alp;
    imod = fn;
}
/*###############################################################################*/ 


#endif