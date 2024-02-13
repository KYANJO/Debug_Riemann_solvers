#define maxiter 1
// #define PI 3.14159265358979323846
// #define deg2rad (PI / 180.0)
// #define rad2deg (180.0 / PI)


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


/* function prototypes */
void riemanntype(double hL, double hR, double uL, double uR, double *hm, double *s1m, double *s2m, bool *rare1, bool *rare2);

void riemann_aug_JCP(int meqn, int mwaves, double hL,
    double hR, double huL, double huR, double hvL, double hvR, 
    double bL, double bR, double uL, double uR, double vL, 
    double vR, double phiL, double phiR, double sE1, double sE2, double* sw, double* fw, int ix, int iy, int idir);

// void riemann_aug_JCP(int meqn, int mwaves, double hL,
//         double hR, double huL, double huR, double hvL, double hvR, 
//         double bL, double bR, double uL, double uR, double vL, 
//         double vR, double phiL, double phiR, double sE1, double sE2, double* sw1, double* sw2, double* sw3, double* fw11, double* fw12, double* fw13, double* fw21, double* fw22, double* fw23, double* fw31, double* fw32, double* fw33, int ix, int iy, int idir);
void cudaflood_rpn2(int idir, int meqn, int mwaves,
                                int maux, double ql[], double qr[],
                                double auxl[], double auxr[],
                                double fwave[], double s[], 
                                double amd_q[], double apd_q[], int ix, int iy,double s_grav, double drytol, double earth_radius, int mcapa);

void cudaflood_rpt2(int idir, int meqn, int mwaves, int maux,
                double q_l[], double q_r[], double aux1[], 
                double aux2[], double aux3[], int imp, 
                double asdq[], double bmasd_q[], double bpasd_q[], int ix, int iy,double s_grav, double drytol, double earth_radius, int mcapa);

/* Main function*/
int main() {
    int meqn = 3, mwaves = 3, maux = 1, mcapa = 0;
    int idir = 0, imp = 0;
    double s_grav = 9.81;
    double drytol = 1.0e-3;
    double earth_radius = 6371000.0;
    int coordinate_system = 1;
    int mx = 16, mbc = 2;
    int ix = 0, iy = 0;
    int my = mx;

    double *ql = (double*)malloc(meqn*sizeof(double));
    double *qr = (double*)malloc(meqn*sizeof(double));
    double *auxl = (double*)malloc(maux*sizeof(double));
    double *auxr = (double*)malloc(maux*sizeof(double));
    double *fwave = (double*)malloc(mwaves*meqn*sizeof(double));
    double *s = (double*)malloc(mwaves*sizeof(double));
    double *amdq = (double*)malloc(mwaves*sizeof(double));
    double *apdq = (double*)malloc(mwaves*sizeof(double));
    
    double *aux1 = (double*)malloc(maux*sizeof(double));
    double *aux2 = (double*)malloc(maux*sizeof(double));
    double *aux3 = (double*)malloc(maux*sizeof(double));
    double *asdq = (double*)malloc(mwaves*sizeof(double));
    double *bmasdq = (double*)malloc(mwaves*sizeof(double));
    double *bpasdq = (double*)malloc(mwaves*sizeof(double));
    
    aux1[0] = 0.0;
    aux2[0] = 0.0;
    aux3[0] = 0.0;
    
    // ql[0] = 0.4473496224909771e2;  
    // ql[1] = -0.1018126050977528e-13;
    // ql[2] = -0.6489769399271391e-9;

    // qr[0] = 0.2519771920562447e2;
    // qr[1] = -0.2016857953898484e-13;
    // qr[2] = -0.3325890167957528e-9;

    // auxl[0] = 0.9999999999996810e2;
    // auxr[0] = 0.9999999999996810e4;

    ql[0] = 2;
    ql[1] = 0.0;
    ql[2] = 0.0;

    qr[0] = 1;
    qr[1] = 0.0;
    qr[2] = 0.0;

    auxl[0] = 0.0;
    auxr[0] = 0.0;
    
    double t0 = 0, Tfinal = 10;
    // double ax = 953236, bx = 959554, ay = 1832407.25, by = 1848572.75;
    double ax = -4, bx = 4, ay = -2, by = 2;
    double dx = (bx- ax)/mx;
    double cfl = 0.75;
    
   double *h = (double*)malloc(mx*sizeof(double));
   double *hu = (double*)malloc(mx*sizeof(double));
   double *hv = (double*)malloc(mx*sizeof(double));
   double *qoldl = (double*)malloc(mwaves*sizeof(double));
   double *qoldr = (double*)malloc(mwaves*sizeof(double));
   double *cxx = (double*)malloc(mx*sizeof(double));
   double *cyy = (double*)malloc(mx*sizeof(double));
   double *Fm = (double*)malloc(mx*sizeof(double));
   double *Fp = (double*)malloc(mx*sizeof(double));

    qoldl = ql;
    qoldr = qr;
    double maxcfl = 0.0;
    double currenttime = 0.0;
    double maxspeed, dt_est,t;
    double dt = 0.01;
    bool cfl_violated = false;

    while (currenttime < Tfinal) {
       h[0] = qoldl[0];
       hu[0] = qoldl[1];
       hv[0] = qoldl[2];
       h[1] = qoldr[0];
       hu[1] = qoldr[1];
       hv[1] = qoldr[2];
       cfl_violated = false; // Reset the flag for each time step
       for(int i=1; i<mx; i++){

           // normal solver    
           cudaflood_rpn2(idir, meqn, mwaves, maux, qoldl, 
                            qoldr, auxl, auxr, fwave, s, amdq, apdq, ix,
                            iy, s_grav, drytol, earth_radius, mcapa);

            // time stepping for the normal solver
            for (int mw = 0; mw < mwaves; mw++) {
                maxspeed = fabs(s[mw]);
                maxcfl = fmax(maxcfl, maxspeed*dt/dx);
                // if (maxcfl > cfl) maxcfl = cfl;
                // dt_est = maxcfl*dx/fmax(maxspeed,1e-6);
            }
            t = currenttime + dt;
            // dt = fmin(dt_est, Tfinal - currenttime);

            if (maxcfl > cfl) {
                dt *= 0.5; // Reduce dt by a factor
                // maxcfl = cfl;
                cfl_violated = true; // Set the flag
                break; // Break out of the spatial loop
            }

            // Transverse solver
            // if(idir ==0 && imp == 0) asdq = amdq;
            // asdq = apdq;
            // cudaflood_rpt2(idir, meqn, mwaves, maux, qoldl, qoldr, aux1, aux2, aux3, imp, 
            //     asdq, bmasdq, bpasdq, ix, iy, s_grav, drytol,earth_radius, 
            //     mcapa);

            printf("s= %.16f,  %.16f, %.16f\n", s[0], s[1], s[2]);
            printf("fwave = %.16f, %.16f,  %.16f\n", fwave[0], fwave[1], fwave[2]);
            printf("fwave = %.16f, %.16f,  %.16f\n", fwave[3], fwave[4], fwave[5]);
            printf("fwave = %.16f,  %.16f, %.16f\n", fwave[6], fwave[7], fwave[8]);
            printf("amdq = %.16f,  %.16f,  %.16f\n", amdq[0], amdq[1], amdq[2]);
            printf("apdq = %.16f,  %.16f,  %.16f\n\n", apdq[0], apdq[1], apdq[2]);

            // printf("\n\n\nasdq = %.16f,  %.16f,  %.16f\n", asdq[0], asdq[1], asdq[2]);
            // printf("bmasdq = %.16f,  %.16f,  %.16f\n", bmasdq[0], bmasdq[1], bmasdq[2]);
            // printf("bpasdq = %.16f,  %.16f,  %.16f\n\n", bpasdq[0], bpasdq[1], bpasdq[2]);
            
            // exit(0);
            // used if we want to use the normal solver only
            h[i] = h[i] - (dt/dx)*(amdq[0] + apdq[0]);
            hu[i] = hu[i] - (dt/dx)*(amdq[1] + apdq[1]);
            hv[i] = hv[i] - (dt/dx)*(amdq[2] + apdq[2]);

            // second order update
            // for (int j = 0; j < mwaves; j++){

            // }
            // Fp = cxx;
            // Fm = cxx;

            // using both normal and transverse solver
            double dy = dx;
            // h[i] = h[i] - (dt/dx)*(amdq[0] + apdq[0]) - (dt/dy)*(bmasdq[0] + bpasdq[0]);
            // hu[i] = hu[i] - (dt/dx)*(amdq[1] + apdq[1]) - (dt/dy)*(bmasdq[1] + bpasdq[1]);
            // hv[i] = hv[i] - (dt/dx)*(amdq[2] + apdq[2]) - (dt/dy)*(bmasdq[2] + bpasdq[2]);
                   // print h[0], h[1]
            // printf("h[0] = %.16f, h[1] = %.16f\n", h[0], h[1]);
            // exit(0);
            // update
            qoldl[0] = h[i-1];
            qoldl[1] = hu[i-1];
            qoldl[2] = hv[i-1];
            
            qoldr[0] = h[i];
            qoldr[1] = hu[i];
            qoldr[2] = hv[i];
            
            printf("dt = %.16f, maxcfl = %.16f, t = %.16f\n", dt, maxcfl,t);
            printf("qoldl = %.16f, %.16f, %.16f\n", qoldl[0],qoldl[1],qoldl[2]);
            printf("qoldr = %.16f, %.16f, %.16f\n\n", qoldr[0],qoldr[1],qoldr[2]);
            // printf("dt = %.16f, dx = %.16f, dt/dx = %.16f, amdq[1] = %.16f, apdq[1] = %.16f\n\n", dt, dx, dt/dx, amdq[1], apdq[1]);
       }
    
       if (cfl_violated) {
            if (t < Tfinal) {
                continue;
            } else {
                printf("CFL violation at time t = %.16f\n", t);
                exit(1);
            }
       }
       
       currenttime += dt;
       dt_est = maxcfl * dx / (maxspeed + 1e-10); // Recalculate dt based on the latest maxspeed
       dt = fmin(dt_est, Tfinal - currenttime);

   }
    
    /* Free allocated memory */
    free(ql); free(qr); free(auxl); free(auxr); free(fwave); free(s); free(amdq); free(apdq);
     free(h); free(hu); free(hv); free(qoldl); free(qoldr); free(aux1); free(aux2); free(aux3); free(asdq); free(bmasdq); free(bpasdq); free(cxx); free(cyy); free(Fm); free(Fp);

    return 0;
} /* end of the main function */

/* Normal Riemann solver for the 2d shallow water equations with topography */
// void cudaflood_rpn2(int idir, int meqn, int mwaves,
//                                 int maux, double ql[], double qr[],
//                                 double auxl[], double auxr[],
//                                 double fwave[], double s[], 
//                                 double amdq[], double apdq[], int ix, int iy)
void cudaflood_rpn2(int idir, int meqn, int mwaves,
                                int maux, double ql[], double qr[],
                                double auxl[], double auxr[],
                                double fwave[], double s[], 
                                double amdq[], double apdq[], int ix, int iy,double s_grav, double drytol, double earth_radius, int mcapa)
{
   
    /* Local variables */
    double wall[3], fw[9], sw[3];
    double hR, hL, huR, huL, hvR, hvL, uR, uL, vR, vL, phiR, phiL;
    double bR, bL, sL, sR, sRoe1, sRoe2, sE1, sE2, uhat, chat;
    double hstar, hstartest, dxdc;
    double s1m, s2m;
    bool rare1, rare2;
    int mw, mu, mv;

    // double pi = 4.0*atan(1.0);
    // double deg2rad = pi/180.0;

    // int N = 4; 
    double deg2rad = 4.0*atan(1.0)/180.0;
      bool debug;
      if (idir == 0 && ix == 7 && iy == 15)
    // if (idir == 0)
      {
        debug = 1;
      }
      else{
        debug = 0;
      }

    /* === Initializing === */
    /* inform of a bad riemann problem from the start */
    if ((qr[0] < 0.0) || (ql[0] < 0.0)) {
        printf("Negative input: hl, hr = %f,%f\n", ql[0], qr[0]);
    }

    /* initialize Riemann problem for grid interface */
    for (mw = 0; mw < mwaves; mw++){
        s[mw] = 0.0;
        fwave[mw] = 0.0;
        fwave[mw + mwaves] = 0.0;
        fwave[mw + 2*mwaves] = 0.0;
    }

    /* set normal direction */
    mu = 1+idir;
    mv = 2-idir;

    /* zero (small) negative values if they exist */
    // left state
    if (qr[0] < 0.0) {
        qr[0] = 0.0;
        qr[1] = 0.0;
        qr[2] = 0.0;
    }

    // right state
    if (ql[0] < 0.0) {
        ql[0] = 0.0;
        ql[1] = 0.0;
        ql[2] = 0.0;
    }

    // if (debug){
    //     printf("ix = %d, iy = %d\n " \ 
    //     "qr[0] = %.16f, ql[0] = %.16f\n" \
    //     "qr[1] = %.16f, ql[1] = %.16f\n" \
    //     "qr[2] = %.16f, ql[2] = %.16f\n\n", ix,iy,qr[0],ql[0],qr[1],ql[1],qr[2],ql[2]);
    // }

    // Skip problem if in a completely dry area
    // if (qr[0] <= drytol && ql[0] <= drytol) {
    //     goto label30;
    // }

    if (ql[0] > drytol || qr[0] > drytol) {
        /* Riemann problem variables */
        // hL = ql[0];
        // hR = qr[0];
        // huL = ql[mu];
        // huR = qr[mu];
        // bL = auxl[0];
        // bR = auxr[0];

        hvL = ql[mv];
        hvR = qr[mv];

        hL = qr[0];
        hR = ql[0];
        huL = qr[mu];
        huR = ql[mu];
        bL = auxr[0];
        bR = auxl[0];

        hvL = qr[mv];
        hvR = ql[mv];

        // Check for wet/dry left boundary
        if (hR > drytol) {
            uR = huR / hR;
            vR = hvR / hR;
            phiR = 0.5 * s_grav * (hR * hR) + (huR * huR) / hR;
        } else {
            hR = 0.0;
            huR = 0.0;
            hvR = 0.0;
            uR = 0.0;
            vR = 0.0;
            phiR = 0.0;
        }

        // Check for wet/dry right boundary
        if (hL > drytol) {
            uL = huL / hL;
            vL = hvL / hL;
            phiL = 0.5 * s_grav * (hL * hL) + (huL * huL) / hL;
        } else {
            hL = 0.0;
            huL = 0.0;
            hvL = 0.0;
            uL = 0.0;
            vL = 0.0;
            phiL = 0.0;
        }
    
        // if (debug){
        //     printf("ix = %d, iy = %d\n " \ 
        //     "hL = %.16f, hR = %.16f\n" \
        //     "huL = %.16f, huR = %.16f\n" \
        //     "hvL = %.16f, hvR = %.16f\n" \
        //     "uL = %.16f, uR = %.16f\n" \
        //     "vL = %.16f, vR = %.16f\n" \
        //     "phiL = %.16f, phiR = %.16f\n" \
        //     "bL = %.16f, bR = %.16f\n\n", ix,iy,hL,hR,huL,huR,hvL,hvR,uL,uR,vL,vR,phiL,phiR,bL,bR);
        // }

        /* left and right surfaces depth inrelation to topography */
        wall[0] = 1.0;
        wall[1] = 1.0;
        wall[2] = 1.0;
        if (hR <= drytol) {
            /* determine the wave structure */
            riemanntype(hL, hL, uL, -uL, &hstar, &s1m, &s2m, &rare1, &rare2);
            // riemann_type(hL, hL, -uL, uL, hstar, s1m, s2m, rare1, rare2);

        //     if (debug){
        //     printf("ix = %d, iy = %d\n " \ 
        //     "hL = %.16f, uL = %.16f\n" \
        //     "hstar = %.16f\n" \
        //     "s1m = %.16f, s2m = %.16f\n" \
        //     "rare1 = %d, rare2 = %d\n\n", ix,iy,hL,uL,hstar,s1m,s2m,rare1,rare2);
        // }

            hstartest = fmax(hL,hstar);
            if (hstartest + bL < bR) {
                /* hL+bL < bR and hstar+bL < bR, so water can't overtop right cell 
                (move into right cell) so right state should become ghost values 
                that mirror left for wall problem) */
                wall[1] = 0.0;
                wall[2] = 0.0;
                hR = hL;
                huR = -huL;
                bR = bL;
                phiR = phiL;
                uR = -uL;
                vR = vL;
                /* here we already have huR =- huL, so we don't need to change it */
            } else if (hL+bL < bR) {
                /* hL+bL < bR and hstar+bL >bR, so we set bR to the water level in 
                the left cell so that water can possibly overtop the right cell (move into the right cell) */ 
                bR = hL + bL;
            }
        } else if (hL <= drytol) { /* right surface is lower than left topo */
            /* determine the Riemann structure */
            riemanntype(hR, hR, -uR, uR, &hstar, &s1m, &s2m, &rare1, &rare2);
            // riemann_type(hR, hR, uR, -uR, hstar, s1m, s2m, rare1, rare2);
            hstartest = fmax(hR,hstar);

            // if (debug){
            //     printf("ix = %d, iy = %d\n " \ 
            //     "hR = %.16f, uR = %.16f\n" \
            //     "hstar = %.16f\n" \
            //     "s1m = %.16f, s2m = %.16f\n" \
            //     "rare1 = %d, rare2 = %d\n\n", ix,iy,hR,uR,hstar,s1m,s2m,rare1,rare2);
            // }

            if (hstartest + bR < bL) //left state should become ghost values that mirror right for wall problem
            {
                wall[0] = 0.0;
                wall[1] = 0.0;
                hL = hR;
                huL = -huR;
                bL = bR;
                phiL = phiR;
                uL = -uR;
                vL = vR;
            } else if (hR+bR < bL) {
                bL = hR + bR;
            }
        }
        // swap variables
        // double htemp = hL;
        // double uhtemp = huL;
        // double vhtemp = hvL;
        // double btemp = bL;
        // double phitemp = phiL;
        // hL = hR;
        // hR = htemp;
        // huL = huR;
        // huR = uhtemp;
        // hvL = hvR;
        // hvR = vhtemp;
        // bL = bR;
        // bR = btemp;
        // phiL = phiR;
        // phiR = phitemp;


        // if (debug){
        //     printf("ix = %d, iy = %d\n " \ 
        //     "hL = %.16f, hR = %.16f\n" \
        //     "huL = %.16f, huR = %.16f\n" \
        //     "hvL = %.16f, hvR = %.16f\n" \
        //     "uL = %.16f, uR = %.16f\n" \
        //     "vL = %.16f, vR = %.16f\n" \
        //     "phiL = %.16f, phiR = %.16f\n" \
        //     "bL = %.16f, bR = %.16f\n\n", ix,iy,hL,hR,huL,huR,hvL,hvR,uL,uR,vL,vR,phiL,phiR,bL,bR);
        // }

        /* determine wave speeds */
        sL = uL - sqrt(s_grav*hL); // 1 wave speed of left state
        sR = uR + sqrt(s_grav*hR); // 2 wave speed of right state

        uhat = (sqrt(s_grav*hL)*uL + sqrt(s_grav*hR)*uR)/(sqrt(s_grav*hL) + sqrt(s_grav*hR)); // Roe average
        chat = sqrt(0.5*s_grav*(hL+hR)); // Roe average
        sRoe1 = uhat - chat; // Roe wave speed 1 wave
        sRoe2 = uhat + chat; // Roe wave speed 2 wave

        sE1 = fmin(sL,sRoe1); // Einfeldt wave speed 1 wave
        sE2 = fmax(sR,sRoe2); // Einfeldt wave speed 2 wave

        // if (ix == 7 && iy == 15) {
        //     if (debug){
        //         printf("ix = %d, iy = %d\n " \ 
        //         "sL = %.16f, sR = %.16f\n" \ 
        //         "sRoe1 = %.16f, sRoe2 = %.16f\n" \ 
        //         "sE1 = %.16f, sE2 = %.16f\n\n", ix,iy,sL,sR,sRoe1,sRoe2,sE1,sE2);
        //     }
        // }
        /* --- end of initializing --- */

        /* === solve Riemann problem === */
        riemann_aug_JCP(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,sw,fw,ix,iy,idir);

        // eliminate ghost fluxes for wall    
        fw[0] *= wall[0];
        fw[mu] *= wall[0];
        fw[mv] *= wall[0];
        sw[0] *= wall[0];

        fw[mwaves + 0] *= wall[1];
        fw[mu + mwaves] *= wall[1];
        fw[mv + mwaves] *= wall[1];
        sw[mu] *= wall[1];

        fw[2*mwaves + 0] *= wall[2];
        fw[2*mwaves + mu] *= wall[2];
        fw[2*mwaves + mv] *= wall[2];
        sw[mv] *= wall[2];

        /* update fwave and corresponding speeds */
        fwave[0] = fw[0];
        fwave[mu] = fw[mu];
        fwave[mv] = fw[mv];
        s[0] = sw[0];

        fwave[mwaves + 0] = fw[mwaves + 0];
        fwave[mu + mwaves] = fw[mu + mwaves];
        fwave[mv + mwaves] = fw[mv + mwaves];
        s[mu] = sw[mu];
       
        fwave[2*mwaves + 0] = fw[2*mwaves + 0];
        fwave[2*mwaves + mu] = fw[2*mwaves + mu];
        fwave[2*mwaves + mv] = fw[2*mwaves + mv];
        s[mv] = sw[mv];
       
    }

    // label30: // (similar to 30 continue in Fortran)

    /* --- Capacity or Mapping from Latitude Longitude to physical space ----*/
    if (mcapa > 0) {
        if (idir == 0) {
            dxdc = earth_radius*deg2rad;
        } else {
            dxdc = earth_radius*cos(auxr[2])*deg2rad;
            // printf("deg2rad =%.16f, dxdc =%.16f\n",deg2rad,dxdc);
        }
        
        // update fwave and corresponding speeds
        for (mw=0; mw<mwaves; mw++) {
            s[mw] *= dxdc;
            fwave[mw] *= dxdc;
            fwave[mw + mwaves] *= dxdc;
            fwave[mw + 2*mwaves] *= dxdc;
        }
    }

     // Debugging
    // if (debug) {
    //     printf("ix = %d, iy = %d\n " \
    //     "s[0] = %.16f, s[1] = %.16f, s[2] = %.16f\n" \
    //     "fwave[0] = %.16f, fwave[1] = %.16f, fwave[2] = %.16f\n" \
    //     "fwave[3] = %.16f, fwave[4] = %.16f, fwave[5] = %.16f\n" \
    //     "fwave[6] = %.16f, fwave[7] = %.16f, fwave[8] = %.16f\n\n", ix,iy,s[0],s[1],s[2],fwave[0],fwave[1],fwave[2],fwave[3],fwave[4],fwave[5],fwave[6],fwave[7],fwave[8]);
    // }

    /* --- compute fluctuations --- */
    amdq[0] = 0.0;
    amdq[1] = 0.0;
    amdq[2] = 0.0;
    apdq[0] = 0.0;
    apdq[1] = 0.0;
    apdq[2] = 0.0;
    int idx; /* mw = idx/3 */
    for (idx = 0; idx < mwaves*3; idx++) {
        if (s[idx/3] < 0.0) { 
            amdq[idx%3] += fwave[idx];  
        } else if (s[idx/3] > 0.0) {
            apdq[idx%3] += fwave[idx]; 
        } else {
            amdq[idx%3] += 0.5 * fwave[idx];  
            apdq[idx%3] += 0.5 * fwave[idx]; 
        }
    }
    // int i;
    // for(mw = 0; mw<mwaves; mw++){
    //     for (i = 0; i<3; i++){
    //         if (s[mw] < 0.0) { 
    //             amdq[i] += fwave[mw*3 + i];  
    //         } else if (s[mw] > 0.0) {
    //             apdq[i] += fwave[mw*3 + i]; 
    //         } else {
    //             amdq[i] += 0.5 * fwave[mw*3 + i];  
    //             apdq[i] += 0.5 * fwave[mw*3 + i]; 
    //         }
    //     }
    // }
}


/* Transverse Riemann solver for the 2d shallow water equations with topography 
@desc: Using The Jacobian matrix from left cell (imp == 0) or right cell (imp == 1) to compute the transverse fluxes.
*/

void cudaflood_rpt2(int idir, int meqn, int mwaves, int maux,
                double ql[], double qr[], double aux1[], 
                double aux2[], double aux3[], int imp, 
                double asdq[], double bmasdq[], double bpasdq[], int ix, int iy,double s_grav, double drytol, double earth_radius, int mcapa) 
{
    /* Access the __constant__ variables in variables.h */
    // double s_grav = 9.81;
    // double drytol = 0.0001;
    // double earth_radius = d_geofloodVars.earth_radius;
    int coordinate_system = 2;
    double pi = 4.0*atan(1.0);
    double deg2rad = pi/180.0;

     int mw, mu, mv;
    double s[3], beta[3];
    double r[3][3];
    double h, u, v;
    double delf1, delf2, delf3;
    double dxdcm, dxdcp, topo1, topo3, eta;

    /* Swapping left to right  (cudaclaw_flux2.cu)*/
    // double *qr = q_l;
    // double *ql = q_r;
    // double *bmasdq = bpasd_q;
    // double *bpasdq = bmasd_q;

    mu = 1+idir;
    mv = 2-idir;

    /* intialize  all components to 0*/
    bmasdq[0] = 0.0;
    bmasdq[mu] = 0.0;
    bmasdq[mv] = 0.0;
    bpasdq[0] = 0.0;
    bpasdq[mu] = 0.0;
    bpasdq[mv] = 0.0;

    h = (imp == 0) ? ql[0] : qr[0];

    bool debug = (idir == 0) ? 1 : 0;
  
    if (h <= drytol) return; // skip problem if dry cell (leaves bmadsq(:) = bpasdq(:) = 0)
    // if (h > drytol) {  
        /* Compute velocities in relevant cell, and other quantities */
        if (imp == 0) {
            // fluctuations being split is left-going
            u = ql[mu] / h;
            v = ql[mv] / h;
            eta = h + aux2[0];
            topo1 = aux1[0];
            topo3 = aux3[0];
        } else {
            // fluctuations being split is right-going
            u = qr[mu] / h;
            v = qr[mv] / h;
            eta = h + aux2[0];
            topo1 = aux1[0];
            topo3 = aux3[0];
        }

        // Debugging check for NaNs
        // if (tid == 0) {
        //     printf("h = %e, u = %e, v = %e, eta = %e, topo1 = %e, topo3 = %e\n", h, u, v, eta, topo1, topo3);
        // }

        /* Check if cell that transverse wave go into are both too high: */
        if (eta < fmin(topo1, topo3)) return; 
        // if (eta >= fmin(topo1, topo3)) {

            /* Check if cell that transverse waves go into are both to high, if so,
            do the splitting (no dry cells), and compute necessary quantities */
            if (coordinate_system == 2) {
                // On the sphere
                if (idir == 1) {
                    dxdcp = earth_radius * deg2rad;
                    dxdcm = dxdcp;
                } else {
                    if (imp == 0) {
                        dxdcp = earth_radius * cos(aux3[2]) * deg2rad;
                        dxdcm = earth_radius * cos(aux1[2]) * deg2rad;
                    } else {
                        dxdcp = earth_radius * cos(aux3[2]) * deg2rad;
                        dxdcm = earth_radius * cos(aux1[2]) * deg2rad;
                    }
                }
            } else {
                // Cartesian
                dxdcp = 1.0;
                dxdcm = 1.0;
            }

            /* Compute some speeds necessary for the Jacobian 
            - Computing upgoing, downgoing waves either in cell on left (if imp==0)
                or on the right (if imp==1) 
            - To achieve this we use q values in cells above and below, however these
                aren't available (only in aux values)
            */
            s[0] = v - sqrt(s_grav * h);
            s[1] = v;
            s[2] = v + sqrt(s_grav * h);

            // Debugging check for NaNs
            // if (ix == 7 && iy == 15) {
            //     // if ((hL >= 0.3280909317849093) && (hR >= 0.3280909317849093)){
            //     if (debug){
            //         printf("ix = %d, iy = %d\n " \
            //         "s[0] = %e, s[1] = %e, s[2] = %e\n", ix,iy, s[0], s[1], s[2]);
            //     }
            // }

            /* Determine asdq decomposition (beta) */
            delf1 = asdq[0];
            delf2 = asdq[mu];
            delf3 = asdq[mv];

            // Debugging check for NaNs
            // if (debug) {
            //     printf("delf1 = %e, delf2 = %e, delf3 = %e\n", delf1, delf2, delf3);
            // }

            beta[0] = (s[2]*delf1 - delf3) / (s[2] - s[0]);
            beta[1] = -u*delf1 + delf2;
            beta[2] = (delf3 - s[0]*delf1) / (s[2] - s[0]);

            /* set-up eigenvectors */
            // r[0] = 1.0;
            // r[mu] = u;
            // r[mv] = s[0];
            r[0][0] = 1.0;
            r[1][0] = u;
            r[2][0] = s[0];

            // r[0 + mwaves] = 0.0;
            // r[mu + mwaves] = 1.0;
            // r[mv + mwaves] = 0.0;
            r[0][1] = 0.0;
            r[1][1] = 1.0;
            r[2][1] = 0.0;

            // r[0 + 2*mwaves] = 1.0;
            // r[mu + 2*mwaves] = u;
            // r[mv + 2*mwaves] = s[2];
            r[0][2] = 1.0;
            r[1][2] = u;
            r[2][2] = s[2];

            // Debugging check for NaNs
            // if (tid == 0) {

            //     printf("beta[0] = %e, beta[1] = %e, beta[2] = %e\n", beta[0], beta[1], beta[2]);
            // }
            // print s, beta, r
            printf("s[0] = %e, s[1] = %e, s[2] = %e\n", s[0], s[1], s[2]);
            printf("beta[0] = %e, beta[1] = %e, beta[2] = %e\n", beta[0], beta[1], beta[2]);
            printf("r[0][0] = %e, r[0][1] = %e, r[0][2] = %e\n", r[0][0], r[0][1], r[0][2]);
            printf("r[1][0] = %e, r[1][1] = %e, r[1][2] = %e\n", r[1][0], r[1][1], r[1][2]);
            printf("r[2][0] = %e, r[2][1] = %e, r[2][2] = %e\n\n", r[2][0], r[2][1], r[2][2]);


            /* Compute transverse fluctuations */
            for (mw = 0; mw < 3; mw++) {
                if ((s[mw] < 0.0) && (eta >= topo1)) {
                    // bmasdq[0] += dxdcm * s[mw]*beta[mw]*r[mw + mwaves];
                    // bmasdq[mu] += dxdcm * s[mw]*beta[mw]*r[mw + mwaves];
                    // bmasdq[mv] += dxdcm * s[mw]*beta[mw]*r[mw + 2*mwaves];
                    bmasdq[0] += dxdcm * s[mw]*beta[mw]*r[0][mw];
                    bmasdq[mu] += dxdcm * s[mw]*beta[mw]*r[1][mw];
                    bmasdq[mv] += dxdcm * s[mw]*beta[mw]*r[2][mw];
                } else if ((s[mw] > 0.0) && (eta >= topo3)) {
                    // bpasdq[0] += dxdcp * s[mw]*beta[mw]*r[mw + mwaves];
                    // bpasdq[mu] += dxdcp * s[mw]*beta[mw]*r[mw + mwaves];
                    // bpasdq[mv] += dxdcp * s[mw]*beta[mw]*r[mw + 2*mwaves];
                    bpasdq[0] += dxdcp * s[mw]*beta[mw]*r[0][mw];
                    bpasdq[mu] += dxdcp * s[mw]*beta[mw]*r[1][mw];
                    bpasdq[mv] += dxdcp * s[mw]*beta[mw]*r[2][mw];
                }
            }
            // Debugging check for NaNs
            // if (tid == 0) {
            //     printf("bmasdq[0] = %e, bmasdq[1] = %e, bmasdq[2] = %e\n", bmasdq[0], bmasdq[1], bmasdq[2]);
            //     printf("bpasdq[0] = %e, bpasdq[1] = %e, bpasdq[2] = %e\n", bpasdq[0], bpasdq[1], bpasdq[2]);
            // }
        // }
    // }
}


/* === Begin fuction riemann_aug_JCP======================================================== @description: - Solves swe give single left and right states
@note: - To use the original solver call with maxiter=1.
       - This solver allows iteration when maxiter > 1. The iteration seems to help  
         with instabilities that arise (with any solver) as flow becomes transcritical 
         over variable topography due to loss of hyperbolicity. 
*/

void riemann_aug_JCP(int meqn, int mwaves, double hL,
    double hR, double huL, double huR, double hvL, double hvR, 
    double bL, double bR, double uL, double uR, double vL, 
    double vR, double phiL, double phiR, double sE1, double sE2, double* sw, double* fw, int ix, int iy, int idir)
// void riemann_aug_JCP(int meqn, int mwaves, double hL,
//         double hR, double huL, double huR, double hvL, double hvR, 
//         double bL, double bR, double uL, double uR, double vL, 
//         double vR, double phiL, double phiR, double sE1, double sE2, double* sw1, double* sw2, double* sw3, double* fw11, double* fw12, double* fw13, double* fw21, double* fw22, double* fw23, double* fw31, double* fw32, double* fw33, int ix, int iy, int idir)
{
    /* Access the __constant__ variables in variables.h */
    double s_grav = 9.81;
    double drytol = 0.0001;

    /* Local variables */
    // double A[9], r[9], lambda[3], del[3], beta[3];
    double lambda[3], beta[3],del[3];
    double A[3][3], r[3][3];
    double delh, delhu, delphi, delb, delnorm;
    double rare1st, rare2st, sdelta, raremin, raremax;
    double criticaltol, convergencetol;
    double criticaltol_2, hustar_interface;
    double s1s2bar, s1s2tilde, hbar, hLstar, hRstar;
    double huRstar, huLstar, uRstar, uLstar, hstarHLL;
    double deldelh, deldelphi;
    double s1m, s2m, hm;
    double det1, det2, det3, determinant;
    bool rare1, rare2, rarecorrector, rarecorrectortest, sonic;
    int mw, k, iter;

    // int mu = 1; // x-direction
    // int mv = 2; // y-direction
    int mu = 1+idir;
    int mv = 2-idir;

    bool debug;
    // if (idir == 0 && (ix-1) == 7 && (iy-1) == 15)
    // if (idir == 0 && ix == 7 && iy == 15)
    if (idir == 0)
    {
      debug = 1;
    }
    else{
      debug = 0;
    }

    /* determine del vectors */
    delh = hR - hL;
    delhu = huR - huL;
    delphi = phiR - phiL;
    delb = bR - bL;
    delnorm = delh * delh + delphi * delphi;

    /* Determine the Riemann structure */
    riemanntype(hL,hR,uL,uR,&hm,&s1m,&s2m,&rare1,&rare2);
    // riemann_type(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2);

    // if (ix == 7 && iy == 0) {
        // if ((hL >= 0.3280909317849093) && (hR >= 0.3280909317849093)){
            // if (debug){
                printf("ix = %d, iy = %d\n " \ 
                "hL = %.16f, hR = %.16f\n" \
                "uL = %.16f, uR = %.16f\n" \
                "hm = %.16f\n" \
                "s1m = %.16f, s2m = %.16f\n" \
                "g = %.16f, drytol = %.16f\n" \
                "rare1 = %d, rare2 = %d\n\n", ix,iy,hL,hR,uL,uR,hm,s1m,s2m,s_grav,drytol,rare1,rare2);
            // }
        // }
    // }
   
    /* For the solver to handle depth negativity, depth dh is included in the decompostion which gives as acess to using the depth positive semidefinite solver (HLLE). This makes the system to have 3 waves instead of 2. where the 1st and 3rd are the eigenpairs are related to the flux Jacobian matrix of the original SWE (since s1<s2<s3, and have been modified by Einfeldt to handle depth non-negativity) and the 2nd is refered to as the the entropy corrector wave since its introduced to correct entropy violating solutions with only 2 waves. */
    
    /* The 1st and 3rd speeds are the eigenvalues of the Jacobian matrix of the original SWE modified by Einfeldt's for use with the HLLE solver. */
    lambda[0] = fmin(sE1, s2m); /* sE1 - flux Jacobian eigen value s2m - Roe speed */
    lambda[2] = fmax(sE2, s1m); /* sE2 - flux Jacobian eigen value s1m - Roe speed*/

    /* Einfeldt's speeds */
    sE1 = lambda[0]; 
    sE2 = lambda[2];

    /* The 2nd speed is the entropy corrector wave speed. */
    lambda[1] = 0.0; /* no strong or significant rarefaction waves */
    
    /* determine the middle state in the HLLE solver */
    hstarHLL = fmax((hL*uL - hR*uR + (sE2 * hR) - (sE1 * hL)) / (sE2 - sE1), 0.0); /* middle state between the two discontinuities (positive semidefinite depth) */

    /* === determine the middle entropy corrector wave === */
    /* rarecorrectortest = .true. provides a more accurate Riemann solution but is more expensive. This is because a nonlinear Riemann solution with  2 nonlinear waves as a linear Riemann solution 3 (or 2 jump discontionuities to approximate 1 smooth nonlinear rarefaction if it's large). When rarecorrectortest = .false. the approximate solution has only 2 jump discontinuities instead of 3, so its less accurate but faster. */
    rarecorrectortest = false;
    rarecorrector = false;
    if (rarecorrectortest) {
        sdelta = lambda[2] - lambda[0];
        raremin = 0.5; /* indicate a large rarefaction wave but not large */
        raremax = 0.9; /* indicate a very large rarefaction wave */
       /* i.e (the total speed difference between the fastest and slowest wave in the Riemann solution = 0.5) */

        if (rare1 && sE1 * s1m < 0.0) raremin = 0.2;
        if (rare2 && sE2 * s2m < 0.0) raremin = 0.2;

        if (rare1 || rare2) {
            /* check which rarefaction is the strongest */
            rare1st = 3.0 * (sqrt(s_grav * hL) - sqrt(s_grav * hm));
            rare2st = 3.0 * (sqrt(s_grav * hR) - sqrt(s_grav * hm));
            if (fmax(rare1st, rare2st) > raremin * sdelta && fmax(rare1st, rare2st) < raremax * sdelta) {
                rarecorrector = true;
                if (rare1st > rare2st) {
                    lambda[1] = s1m;
                } else if (rare2st > rare1st) {
                    lambda[1] = s2m;
                } else {
                    lambda[1] = 0.5 * (s1m + s2m);
                }
            }
        }
        if (hstarHLL < fmin(hL, hR) / 5.0) rarecorrector = false;
    }

    /* determining modified eigen vectors */
    for (mw = 0; mw < mwaves; mw++) {   
        // r[mw] = 1.0; 
        // r[mw + mwaves] = lambda[mw]; 
        // r[mw + 2*mwaves] = pow(lambda[mw],2.0);
        r[0][mw] = 1.0;
        r[1][mw] = lambda[mw];
        r[2][mw] = pow(lambda[mw],2.0);
    }

    /* no strong rarefaction wave */
    if (!rarecorrector) {
        lambda[1]= 0.5*(lambda[0] + lambda[2]);
        // r[mwaves] = 0.0; // r[0,1]
        // r[mwaves + mu] = 0.0; // r[1,1]
        // r[mwaves + mv] = 1.0; // r[2,1]
        r[0][1] = 0.0;
        r[1][1] = 0.0;
        r[2][1] = 1.0;
    }

    /* === Determine the steady state wave === */
    criticaltol = fmax(drytol*s_grav, 1.0e-6);
    criticaltol_2 = sqrt(criticaltol);
    deldelh = -delb;
    deldelphi = -0.5 * (hR + hL) * (s_grav * delb); /* some approximation of the source term \int_{x_{l}}^{x_{r}} -g h b_x dx */

    /* determine a few quantities needed for steady state wave if iterated */
    hLstar = hL;
    hRstar = hR;
    uLstar = uL;
    uRstar = uR;
    huLstar = uLstar * hLstar;
    huRstar = uRstar * hRstar;

    /* iterate to better find the steady state wave */
    convergencetol = 1e-6;
    for (iter=1; iter <= maxiter; iter++) {
        /* determine steady state wave (this will be subtracted from the delta vectors */
        if (fmin(hLstar,hRstar) < drytol && rarecorrector) {
            rarecorrector = false;
            hLstar = hL;
            hRstar = hR;
            uLstar = uL;
            uRstar = uR;
            huLstar = uLstar*hLstar;
            huRstar = uRstar*hRstar;
            lambda[1] = 0.5*(lambda[0] + lambda[2]);
            // r[mwaves] = 0.0; // r[0,1]
            // r[mwaves + mu] = 0.0; // r[1,1]
            // r[mwaves + mv] = 1.0; // r[2,1]
            r[0][1] = 0.0;
            r[1][1] = 0.0;
            r[2][1] = 1.0;
        }

        /* For any two states; Q_i and Q_i-1, eigen values of SWE must satify: lambda(q_i)*lambda(q_i-1) = u^2 -gh, writing this conditon as a function of Q_i and Q_i-1, u and h become averages in lambda(q_i)*lambda(q_i-1) = u^2 -gh and these averages are denoted by bar and tilde. */
        hbar = fmax(0.5 * (hLstar + hRstar), 0.0);
        s1s2bar = 0.25 * pow((uLstar + uRstar),2) - (s_grav * hbar);
        s1s2tilde = fmax(0.0, uLstar * uRstar) - (s_grav * hbar);

        /* Based on the above conditon, smooth staedy state over slopping bathymetry cannot have a sonic point. Therefore, for regions with monotonically varying bathymetry, steady-state flow is either entirely subsonic (-u^2 +gh > 0) or entirely supersonic. */
        sonic = false;
        if (fabs(s1s2bar) <= criticaltol) {
            sonic = true;
        } else if (s1s2bar * s1s2tilde <= criticaltol * criticaltol) {
            sonic = true;
        } else if (s1s2bar * sE1 * sE2 <= criticaltol * criticaltol) {
            sonic = true;
        } else if (fmin(fabs(sE1), fabs(sE2)) < criticaltol_2) {
            sonic = true;
        } else if (sE1 < criticaltol_2 && s1m > -criticaltol_2) {
            sonic = true;
        } else if (sE2 > -criticaltol_2 && s2m < criticaltol_2) {
            sonic = true;
        } else if ((uL + sqrt(s_grav * hL)) * (uR + sqrt(s_grav * hR)) < 0.0) {
            sonic = true;
        } else if ((uL - sqrt(s_grav * hL)) * (uR - sqrt(s_grav * hR)) < 0.0) {
            sonic = true;
        }

        /* find jump in h, deldelh */
        if (sonic) {
            deldelh = -delb;
        } else {
            deldelh = delb * s_grav * hbar / s1s2bar;
        }

        /* find bounds in case of critical state resonance, or negative states */
        if (sE1 < -criticaltol && sE2 > criticaltol) {
            deldelh = fmin(deldelh, hstarHLL * (sE2 - sE1) / sE2);
            deldelh = fmax(deldelh, hstarHLL * (sE2 - sE1) / sE1);
        } else if (sE1 >= criticaltol) {
            deldelh = fmin(deldelh, hstarHLL * (sE2 - sE1) / sE1);
            deldelh = fmax(deldelh, -hL);
        } else if (sE2 <= -criticaltol) {
            deldelh = fmin(deldelh, hR);
            deldelh = fmax(deldelh, hstarHLL * (sE2 - sE1) / sE2);
        }

        /* find jump in phi, ddphi */
        if (sonic) {
            deldelphi = -s_grav * hbar * delb;
        } else {
            deldelphi = -delb * s_grav * hbar * s1s2tilde / s1s2bar;
        }

        /* find bounds in case of critical state resonance, or negative states */
        deldelphi = fmin(deldelphi, s_grav * fmax(-hLstar * delb, -hRstar * delb));
        deldelphi = fmax(deldelphi, s_grav * fmin(-hLstar * delb, -hRstar * delb));

        /* determine the delta vectors */
        del[0] = delh - deldelh;
        del[1] = delhu;
        del[2] = delphi - deldelphi;  

        /* Determine coefficients beta(k) using crammer's rule
          first determine the determinant of the eigenvector matrix */
        // det1 = r[0]*(r[mwaves + mu]*r[2*mwaves + mv] - r[2*mwaves + mu]*r[mwaves + mv]);
        // det2 = r[mwaves]*(r[mu]*r[2*mwaves + mv] - r[2*mwaves + mu]*r[mv]);
        // det3 = r[2*mwaves]*(r[mu]*r[mwaves + mv] - r[mwaves + mu]*r[mv]);
        det1 = r[0][0]*(r[1][1]*r[2][2] - r[1][2]*r[2][1]);
        det2 = r[0][1]*(r[1][0]*r[2][2] - r[1][2]*r[2][0]);
        det3 = r[0][2]*(r[1][0]*r[2][1] - r[1][1]*r[2][0]);
        determinant = det1 - det2 + det3;

        /* solve for beta(k) */
        for(k=0; k < 3; k++)
        {   
            for(mw=0; mw < 3; mw++)
            {
                // A[mw] = r[mw]; 
                // A[mw + mwaves] = r[mw + mwaves];
                // A[mw + 2*mwaves] = r[mw + 2*mwaves];
                A[0][mw] = r[0][mw];
                A[1][mw] = r[1][mw];
                A[2][mw] = r[2][mw];
            }
            // A[k] = del[0];
            // A[mwaves + k] = del[1];
            // A[2*mwaves + k] = del[2];
            // det1 = A[0]*(A[mwaves + mu]*A[2*mwaves + mv] - A[2*mwaves + mu]*A[mwaves + mv]);
            // det2 = A[mwaves]*(A[mu]*A[2*mwaves + mv] - A[2*mwaves + mu]*A[mv]);
            // det3 = A[2*mwaves]*(A[mu]*A[mwaves + mv] - A[mwaves + mu]*A[mv]);
            A[0][k] = del[0];
            A[1][k] = del[1];
            A[2][k] = del[2];
            det1 = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]);
            det2 = A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
            det3 = A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
            beta[k] = (det1 - det2 + det3)/determinant;
        }

        /* exit if things aren't changing */
        if (fabs(pow(del[0],2)+pow(del[2],2.0) - delnorm) < convergencetol) break;

        delnorm = pow(del[0],2)+pow(del[2],2.0); /* update delnorm */

        /* find new states qLstar and qRstar on either side of interface */
        hLstar = hL;
        hRstar = hR;
        uLstar = uL;
        uRstar = uR;
        huLstar = uLstar*hLstar;
        huRstar = uRstar*hRstar;

        /* left state depth and momentum updates */
        for (mw=0; mw < mwaves; mw++)
        {
            if (lambda[mw] < 0.0)
            {
                // hLstar = hLstar + beta[mw]*r[mw]; 
                // huLstar = huLstar + beta[mw]*r[mw + mwaves]; 
               hLstar = hLstar + beta[mw]*r[0][mw];
               huLstar = huLstar + beta[mw]*r[1][mw];
            }
        }

        /* right state depth and momentum updates */
        for (mw = mwaves-1; mw >= 0; mw--)
        {
            if (lambda[mw] > 0.0)
            {
                // hRstar = hRstar - beta[mw]*r[mw]; 
                // huRstar = huRstar - beta[mw]*r[mw + mwaves]; 
                hRstar = hRstar - beta[mw]*r[0][mw];
                huRstar = huRstar - beta[mw]*r[1][mw];
            }
        }

        /* left state velocity update */
        if (hLstar > drytol) 
        {
            uLstar = huLstar/hLstar;
        }
        else  /* dry state */
        {
            hLstar = fmax(hLstar,0.0);
            uLstar = 0.0;
        }

        /* right state velocity update */
        if (hRstar > drytol) 
        {
            uRstar = huRstar/hRstar;
        }
        else /* dry state */
        {
            hRstar = fmax(hRstar,0.0);
            uRstar = 0.0;
        }
    } /* end of  iteration on the Riemann problem*/

    /* === determine the fwaves and speeds=== */
    // for (mw=0; mw < mwaves; mw++)
    // {
    //     sw[mw] = lambda[mw];
    //     // fw[mw] = beta[mw]*r[mw + mwaves]; 
    //     // fw[mw + mwaves] = beta[mw]*r[mw + 2*mwaves]; 
    //     // fw[mw + 2*mwaves] = beta[mw]*r[mw + mwaves]; 
    //     fw[mw] = beta[mw]*r[1][mw];
    //     fw[mw + mwaves] = beta[mw]*r[2][mw];
    //     fw[mw + 2*mwaves] = beta[mw]*r[1][mw];
    // }

    fw[0] = beta[0]*r[1][0];
    fw[mu] = beta[0]*r[2][0];
    fw[mv] = beta[0]*r[1][0];
    sw[0] = lambda[0];

    fw[mwaves + 0] = beta[1]*r[1][1];
    fw[mwaves + mu] = beta[1]*r[2][1];
    fw[mwaves + mv] = beta[1]*r[1][1];
    sw[1] = lambda[1];

    fw[2*mwaves + 0] = beta[2]*r[1][2];
    fw[2*mwaves + mu] = beta[2]*r[2][2];
    fw[2*mwaves + mv] = beta[2]*r[1][2];
    sw[2] = lambda[2];

    // *sw1 = lambda[0];
    // *sw2 = lambda[1];
    // *sw3 = lambda[2];
    // // mw = 0;
    // *fw11 = beta[0]*r[1][0];
    // *fw21 = beta[0]*r[2][0];
    // *fw31 = beta[0]*r[1][0];
    // // mw = 1;
    // *fw12 = beta[1]*r[1][1];
    // *fw22 = beta[1]*r[2][1];
    // *fw32 = beta[1]*r[1][1];
    // // mw = 2;
    // *fw13 = beta[2]*r[1][2];
    // *fw23 = beta[2]*r[2][2];
    // *fw33 = beta[2]*r[1][2];

    // find transverse components (ie huv jumps)
    fw[mv] *= vL;
    fw[2*mwaves + mv] *= vR;
    fw[mwaves + mv] = 0.0;
    // *fw31 *= vL;
    // *fw33 *= vR;
    // *fw32 = 0.0;

    hustar_interface = hL*uL + fw[0];
    if (hustar_interface <= 0.0) {
        fw[mv] += (hR * uR * vR - hL * uL * vL - fw[mv] - fw[2*mwaves + mv]);
    } else {
        fw[2*mwaves + mv] += (hR * uR * vR - hL * uL * vL - fw[mv] - fw[2*mwaves + mv]);
    }
    // hustar_interface = hL*uL + *fw11;
    // if (hustar_interface <= 0.0) {
    //     *fw31 += (hR * uR * vR - hL * uL * vL - *fw31 - *fw33);
    // } else {
    //     *fw33 += (hR * uR * vR - hL * uL * vL - *fw31 - *fw33);
    // }

}

/* === Begin fuction Riemann type ============
 @description: Determines the Riemann structure (wave-type in each family)
*/

void riemanntype(double hL, double hR, double uL, double uR, double *hm, 
                            double *s1m, double *s2m, bool *rare1, bool *rare2)
{
    /* Access the __constant__ variables in variables.h */
    double s_grav = 9.81;
    double drytol = 0.0001;

    // Local variables
    double u1m, u2m, h0, F_max, F_min, dfdh, F0, slope, gL, gR;
    double sqrtgh1, sqrtgh2;
    // double um;
    int iter; 

    // Test for Riemann structure
    double h_min = fmin(hR,hL);
    double h_max = fmax(hR,hL);
    double delu = uR - uL;

    /* Have dry state on either side 
    - Only one rarefaction wave
    - another shock wave has 0 jump and moves at the same speed as one edge of the    rarefaction  wave */
    if (h_min <= drytol)
    {
        *hm = 0.0;
        // um = 0.0;
    
        /* Either hR or hL is almost zero, so the expression below corresponds
           to either Eqn. (54a) or Eqn. (54b) in the JCP paper */
        *s1m = uR + uL - 2.0 * sqrt(s_grav * hR) + 2.0 * sqrt(s_grav * hL);
        *s2m = *s1m; 
        *rare1 = (hL <= 0.0) ? false : true;
        *rare2 = !(*rare1);
    } else {
        F_min = delu + 2.0 * (sqrt(s_grav * h_min) - sqrt(s_grav * h_max));
        F_max = delu + (h_max - h_min) * sqrt(0.5 * s_grav * (h_max + h_min) / (h_max * h_min));

        if (F_min > 0.0){  // 2-rarefactions
            /* Eqn (13.56) in the FVMHP book */
            *hm = (1.0 / (16.0 * s_grav)) * pow(fmax(0.0, -delu + 2.0 * (sqrt(s_grav * hL) + sqrt(s_grav * hR))), 2);
            // um = copysign(1.0, *hm) * (uL + 2.0 * (sqrt(s_grav * hL) - sqrt(s_grav * *hm)));
            *s1m = uL + 2.0 * sqrt(s_grav * hL) - 3.0 * sqrt(s_grav * *hm);
            *s2m = uR - 2.0 * sqrt(s_grav * hR) + 3.0 * sqrt(s_grav * *hm);
            *rare1 = true;
            *rare2 = true;
        } else if (F_max <= 0.0) { // 2-shocks
            /* Below it solves for the intersection of two Hugoniot loci to get the
            accurate Riemann solution */
            /* Root finding using a Newton iteration on sqrt(h) */
            h0 = h_max;
            for (iter = 1; iter <= maxiter; iter++) {
                gL = sqrt(0.5 * s_grav * (1.0 / h0 + 1.0 / hL));
                gR = sqrt(0.5 * s_grav * (1.0 / h0 + 1.0 / hR));
                F0 = delu + (h0 - hL) * gL + (h0 - hR) * gR;
                dfdh = gL - s_grav * (h0 - hL) / (4.0 * h0 * h0 * gL) + gR - s_grav * (h0 - hR) / (4.0 * h0 * h0 * gR);
                slope = 2.0 * sqrt(h0) * dfdh;
                h0 = pow(sqrt(h0) - F0 / slope, 2);
            }
            *hm = h0;
            /* u1m and u2m are Eqns (13.19) and (13.20) in the FVMHP book */
            u1m = uL - (*hm - hL) * sqrt(0.5 * s_grav * (1.0 / *hm + 1.0 / hL));
            u2m = uR + (*hm - hR) * sqrt(0.5 * s_grav * (1.0 / *hm + 1.0 / hR));
            // um = 0.5 * (u1m + u2m);
            *s1m = u1m - sqrt(s_grav * *hm);
            *s2m = u2m + sqrt(s_grav * *hm);
            *rare1 = false;
            *rare2 = false;
        } else { // 1-shock or 1-rarefaction
            h0 = h_min;
            for (iter = 1; iter <= maxiter; iter++) {
                F0 = delu + 2.0 * (sqrt(s_grav * h0) - sqrt(s_grav * h_max)) + (h0 - h_min) * sqrt(0.5 * s_grav * (1.0 / h0 + 1.0 / h_min));
                slope = (F_max - F0) / (h_max - h_min);
                h0 = h0 - F0 / slope;
            }
            *hm = h0;
            sqrtgh2 = sqrt(s_grav * *hm);
            if (hL > hR) {
                sqrtgh1 = sqrt(s_grav * hL);
                /* Eqn (13.55) in the FVMHP book */
                // um = uL + 2.0 * sqrtgh1 - 2.0 * sqrtgh2;
                *s1m = uL + 2.0 * sqrtgh1 - 3.0 * sqrtgh2;
                *s2m = uL + 2.0 * sqrtgh1 - sqrtgh2;

                *rare1 = true;
                *rare2 = false;
            } else {
                sqrtgh1 = sqrt(s_grav * hR);
                // um = uR - 2.0 * sqrtgh1 + 2.0 * sqrtgh2;
                *s1m = uR - 2.0 * sqrtgh1 + sqrtgh2;
                *s2m = uR - 2.0 * sqrtgh1 + 3.0 * sqrtgh2;
                *rare1 = false;
                *rare2 = true;
            }
        }
    }
} /* End of riemanntype function */



//  void riemann_type(double hL, double hR, double uL, double uR, double hm, 
//     double s1m, double s2m, bool rare1, bool rare2)
// {
//     double g = 9.81;
//     double drytol = 0.0001;

//     double um,u1m,u2m,delu;
//     double h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR;
//     double sqrtgh1,sqrtgh2;
//     int iter;

//     /* Test for Riemann structure */
//     h_min = min(hR,hL);
//     h_max = max(hR,hL);
//     delu = uR - uL;

//    if (h_min <= drytol){
//     hm = 0.0;
//     um = 0.0;
//     s1m = uR + uL - 2.0*sqrt(g*hR) + 2.0*sqrt(g*hL);
//     s2m = uR + uL - 2.0*sqrt(g*hR) + 2.0*sqrt(g*hL);

//     if (hL <= 0.0){
//         rare1 = false;
//         rare2 = true;
//     } else {
//         rare1 = true;
//         rare2 = false;
//     }
//    } else {
//         F_min = delu + 2.0*(sqrt(g*h_min) - sqrt(g*h_max));
//         F_max = delu + (h_max - h_min)*(sqrt(0.5*g*(h_max + h_min)/(h_max*h_min)));

//         if (F_min < 0.0) {
//             /* 2-rarefactions */
//             hm = (1.0/(16.0*g))*(pow(max(0.0,-delu + 2.0*(sqrt(g*hL) + sqrt(g*hR))),2));
//             um = copysign(1.0,hm)*(uL + 2.0*(sqrt(g*hL) - sqrt(g*hm)));
//             s1m = uL + 2.0*sqrt(g*hL) - 3.0*sqrt(g*hm);
//             s2m = uR - 2.0*sqrt(g*hR) + 3.0*sqrt(g*hm);
//             rare1 = true;
//             rare2 = true;
//         } else if (F_max <= 0.0) {
//             /* 2-shocks */
//             /* Root finding using a Newton iteration on sqrt(h) */
//             h0 = h_max;
//             for (iter = 1; iter <= maxiter; iter++) {
//                 gL = sqrt(0.5*g*(1.0/h0 + 1.0/hL));
//                 gR = sqrt(0.5*g*(1.0/h0 + 1.0/hR));
//                 F0 = delu + (h0 - hL)*gL + (h0 - hR)*gR;
//                 dfdh = gL - g*(h0 - hL)/(4.0*h0*h0*gL) + gR - g*(h0 - hR)/(4.0*h0*h0*gR);
//                 slope = 2.0*sqrt(h0)*dfdh;
//                 h0 = pow(sqrt(h0) - F0/slope,2);
//             }
//             hm = h0;
//             /* u1m and u2m are Eqns (13.19) and (13.20) in the FVMHP book */
//             u1m = uL - (hm - hL)*sqrt(0.5*g*(1.0/hm + 1.0/hL));
//             u2m = uR + (hm - hR)*sqrt(0.5*g*(1.0/hm + 1.0/hR));
//             um = 0.5*(u1m + u2m);
//             s1m = u1m - sqrt(g*hm);
//             s2m = u2m + sqrt(g*hm);
//             rare1 = false;
//             rare2 = false;
//         } else {
//             /* 1-shock or 1-rarefaction */
//             h0 = h_min;
//             for (iter = 1; iter <= maxiter; iter++) {
//                 F0 = delu + 2.0*(sqrt(g*h0) - sqrt(g*h_max)) + (h0 - h_min)*sqrt(0.5*g*(1.0/h0 + 1.0/h_min));
//                 slope = (F_max - F0)/(h_max - h_min);
//                 h0 = h0 - F0/slope;
//             }
//             hm = h0;
//             sqrtgh2 = sqrt(g*hm);
//             if (hL > hR) {
//                 sqrtgh1 = sqrt(g*hL);
//                 /* Eqn (13.55) in the FVMHP book */
//                 um = uL + 2.0*sqrtgh1 - 2.0*sqrtgh2;
//                 s1m = uL + 2.0*sqrtgh1 - 3.0*sqrtgh2;
//                 s2m = uL + 2.0*sqrtgh1 - sqrtgh2;

//                 rare1 = true;
//                 rare2 = false;
//             } else {
//                 sqrtgh1 = sqrt(g*hR);
//                 um = uR - 2.0*sqrtgh1 + 2.0*sqrtgh2;
//                 s1m = uR - 2.0*sqrtgh1 + sqrtgh2;
//                 s2m = uR - 2.0*sqrtgh1 + 3.0*sqrtgh2;
//                 rare1 = false;
//                 rare2 = true;
//             }
//         }
//    }

// }

