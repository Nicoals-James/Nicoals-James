
#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <math.h>
#define XILLVER11 "xillverD-4.fits"

double abun0,xi0,logden0,cos0,gamma0;//give out specify parameaters to get rid of them 
double gam0[128]={0};
int i,j,k,m,n,it;//index for cycle
#define tbin  128
double flux_fin[tbin][3000];
/*****************you should give out the energy bin(nener), and gamma***********************/

int main(int argc, char*argv[])
{
       
       long ixi,ilogden,icose,igam,iabun,ie,irow,je,quit;//index for cycle
       abun0=1;xi0=2;logden0=17;cos0=45;gamma0=2;
       /*
       ---parameter--------Min-----------Max--
       Abundance_Fe        0.5           20
       logXi                0          4.69897 
       log_Density         15            19
       Inci_cose          18.195       87.124   (degree)
       */


       fitsfile *fptr;
       int status=0,hdutype=2;
       long  frow = 1, felem = 1, nelems, nrow;//locate the data in the table 
       int colnum ;//table column to read  
       int anynul;
       float    float_nulval = 0.;
       static long  ngam=12,nabun=5,nlogden=9,ncose=10,nxi=15;//number of values
       static long ne_loc;
       static float *gam,*abun,*cose,*logxi,*logden;//dynamic array of parameters' value (cose in degree)
       static float *energy0;//original energy array
       static float *emission;
       int nspec;//total number of values
       static double *flux0;//the spectrum in the table
       static double *flux1;//the spectrum in bins you want
       static long nener;//the number of energy bins you want
       nspec = ngam*nabun*nlogden*ncose*nxi;


       int imin,imax;//original index for interpolation
       int iabun0,ixi0,icos0,ilogden0,igam0;//the index of specify parameter
       double abun_tmp1,xi_tmp1,cos_tmp1,den_tmp1,gam_tmp1;//the precise index |___tmp1____*para0*____tmp0_____|
       double abun_tmp0,xi_tmp0,cos_tmp0,den_tmp0,gam_tmp0;//the precise index |___tmp1____*para0*____tmp0_____|
       double abun_tmp[2],xi_tmp[2],cos_tmp[2],den_tmp[2],gam_tmp[2];
       double elow,ehigh;//the range of energy 
       //double y0000,y0001,y0010,y0011,y0100,y0101,y0110,y0111,y1000,y1001,y1010,y1011,y1100,y1101,y1110,y1111;//the interpolation values
       double y[2][2][2][2]={0};//2-D 4-orders tensor for interpolation  

       static double *energy1,*energy2;//new energy bin after rebinning
       double egap_tmp,etmp,etmp1;
       ffopen(&fptr, XILLVER11, READONLY, &status);
/*******************move to the HDU "PARAMETERS"*********************/
       //ffmrhd(fptr, 1, &hdutype, &status);
       //locate_hdu(fptr);
       /*
       row1: Gamma
       row2: Abundance_Fe
       row3: logXi
       row4: Density
       row5: incident angle(cos)
       */
      /*
      ffmrhd(fptr, 1, &hdutype, &status);
      int chdunum=0;
    ffghdn(fptr,&chdunum);
    printf("The current HDU is HDU  %d \n",chdunum);
    locate_hdu(fptr);
        nelems = 1;
        colnum=10;
    printf("%d\n",colnum);
        frow = 1;
       ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &ngam,
        &anynul, &status);
    printf("%ld",ngam);
        frow = 2;
       ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &nabun,
        &anynul, &status);

        frow = 3;
       ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &nxi,
        &anynul, &status);

        frow = 4;
       ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &nlogden,
            &anynul, &status);
     
        frow = 5;
       ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &ncose,
            &anynul, &status);
        *///we can write down the NumbVals directly

//********** free all the dynamic array
        if(gam != NULL){ free((void *) gam); gam = NULL;}
        if(abun != NULL){ free((void *) abun); abun = NULL;}
        if(logxi != NULL){ free((void *) logxi); logxi = NULL;}
        if(logden != NULL){ free((void *) logden); logden = NULL;}
        if(cose != NULL){ free((void *) cose); cose = NULL;}
        if(emission != NULL){ free((void *) cose); cose = NULL;}
        


/********************************read the values of the parameters************************************/
        
        ffmrhd(fptr, 1, &hdutype, &status);
        locate_hdu(fptr);

    //allocate memory to dynamic arrays
        gam = (float *) malloc( ngam * sizeof(float));
        abun = (float *) malloc( nabun * sizeof(float));
        cose = (float *) malloc( ncose * sizeof(float));
        logxi = (float *) malloc( nxi * sizeof(float));
        logden = (float *) malloc( nlogden * sizeof(float));
    //read the value
        
        colnum = 10;
        frow = 1;
        nelems = ngam;
        ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, gam,&anynul, &status);
        
        frow = 2;
        nelems = nabun;
        ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, abun,&anynul, &status);

        frow = 3;  
        nelems = nxi;
        ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, logxi,&anynul, &status);

        frow = 4;
        nelems = nlogden;
        ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, logden,&anynul, &status);

        frow = 5;
        nelems = ncose;
        ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, cose,&anynul, &status);
        
        /* for(i=0;i<ncose;i++){
            cose[i]=cos(cose[i]/180.*M_PI);
        } */
        

       printf("max cos %f \n", cose[0]);


/*******************move to the HDU "ENERGIES", read energy value*********************/

       ffmrhd(fptr, 1, &hdutype, &status);
       locate_hdu(fptr);

       ffgnrw(fptr, &ne_loc, &status);//read energy bins

       ne_loc++;
   
       if(energy0 != NULL){ free((void *) energy0); energy0 = NULL;}
       energy0 = (float *) malloc(ne_loc * sizeof(float));//dynamic array 
       nelems = ne_loc - 1;
       colnum = 1;
       frow = 1;
       ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, energy0,&anynul, &status);
       //printf("3  %ld",energy0);
       nelems = 1;
       colnum = 2;
       frow = ne_loc - 1;
       ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, &energy0[ne_loc - 1], &anynul, &status);
       printf("min energy : %f \n ",energy0[0]);
       printf("max energy : %f \n ",energy0[ne_loc-1]);

/**/
/*****************************move to the HDU "SPECTRA" read the emission************************************/
        long nelements;
        ffmrhd(fptr, 1, &hdutype, &status);
        locate_hdu(fptr);
        
       
        ffgrsz(fptr, &nrow, &status);//read rows 
        nelements = nrow * (ne_loc - 1);
        printf("the optimal number of row to read : %ld \n",nrow);

        emission = (float *) malloc(nspec * (ne_loc - 1) * sizeof(float));//dynamic array
       /*  double et;
        et = (float) malloc(2999*sizeof(float)); */

        colnum = 2;
        
       for (irow = 0; irow < nspec; irow += nrow) {
// the last block to read may be smaller:
    
        if ((nspec - irow) < nrow) nelements = (nspec - irow) * (ne_loc - 1);
        igam = irow / (nabun * nxi * nlogden * ncose) + 1;
        iabun = (irow - (igam - 1) * nabun * nxi * nlogden * ncose) /
                (nxi * nlogden * ncose) + 1;
        ixi = (irow - (igam - 1) * nabun * nxi * nlogden * ncose - 
              (iabun - 1)* nxi * nlogden * ncose) /
              (nlogden * ncose) + 1;
        ilogden = (irow - (igam - 1) * nabun * nxi * nlogden * ncose - 
                (iabun - 1)* nxi * nlogden * ncose - 
                (ixi - 1) * nlogden * ncose) / ncose + 1;
        icose = irow - (igam - 1) * nabun * nxi * nlogden * ncose - 
                      (iabun - 1)* nxi * nlogden * ncose - 
                      (ixi - 1) * nlogden * ncose -
                      (ilogden - 1) * ncose + 1;
                      
       
        
        ffgcv(fptr, TFLOAT, 2, irow + 1, 1, nelements, &float_nulval, 
              &emission[(ne_loc - 1) * (icose - 1) + 
              (ne_loc - 1) * ncose * (ilogden - 1) + 
              (ne_loc - 1) * ncose * nlogden * (ixi - 1) + 
              (ne_loc - 1) * ncose * nlogden * nxi * (iabun - 1) + 
              (ne_loc - 1) * ncose * nlogden * nxi * nabun * (igam - 1)], 
              &anynul, &status);//the sequence of emission is (gamma, abundence, logxi, density, cose)
      }
    printf("emission %f \n",emission[0]);
    ffclos(fptr, &status);
    printf("end of all HDU \n");
/***************************** end of reading xillver fits file ********************************/
        
/*************************interpolate the emission********************************/
//get rid of abun, logxi, logden, cose. and find the index for given gamma
    // given abundance, find the corresponding index in abun[]:
        imin = 0;
        imax = nabun;
        iabun0 = nabun / 2;
        while ((imax - imin) > 1) {
            if (abun0 >= abun[iabun0 - 1]) imin = iabun0;
            else imax = iabun0;
            iabun0 = (imin + imax) / 2;
        }
        if (iabun0 == 0) iabun0 = 1;
        abun_tmp1 = (abun0 - abun[iabun0 - 1]) / (abun[iabun0] - abun[iabun0 - 1]);
        if (abun_tmp1 < 0.) abun_tmp1 = 0.;
        if (abun_tmp1 > 1.) abun_tmp1 = 1.;
        abun_tmp0 = 1. - abun_tmp1;
        abun_tmp[0]=abun_tmp0;
        abun_tmp[1]=abun_tmp1;

    // given logxi, find the corresponding index in logxi[]:
        imin = 0;
        imax = nxi;
        ixi0 = nxi / 2;
        while ((imax - imin) > 1) {
            if (xi0 >= logxi[ixi0 - 1]) imin = ixi0;
            else imax = ixi0;
            ixi0 = (imin + imax) / 2;
        }
        if (ixi0 == 0) ixi0 = 1;
        xi_tmp1 = (xi0 - logxi[ixi0 - 1]) / (logxi[ixi0] - logxi[ixi0 - 1]);
        if (xi_tmp1 < 0.) xi_tmp1 = 0.;
        if (xi_tmp1 > 1.) xi_tmp1 = 1.;
        xi_tmp0 = 1. - xi_tmp1;   
        xi_tmp[0]=xi_tmp0;
        xi_tmp[1]=xi_tmp1;

    // given density, find the corresponding index in logden[]:
        imin = 0;
        imax = nlogden;
        ilogden0 = nlogden / 2;
        while ((imax - imin) > 1) {
            if (logden0 >= logden[ilogden0 - 1]) imin = ilogden0;
            else imax = ilogden0;
            ilogden0 = (imin + imax) / 2;
        }
        if (ilogden0 == 0) ilogden0 = 1;
        den_tmp1 = (logden0 - logden[ilogden0 - 1]) / (logden[ilogden0] - logden[ilogden0 - 1]);
        if (den_tmp1 < 0.) den_tmp1 = 0.;
        if (den_tmp1 > 1.) den_tmp1 = 1.;
        den_tmp0 = 1. - den_tmp1;  
        den_tmp[0]=den_tmp0;
        den_tmp[1]=den_tmp1;

    // given cose, find the corresponding index in cose[]:
        imin = 0;
        imax = ncose;
        icos0 = ncose / 2;
        while ((imax - imin) > 1) {
            if (cos0 >= cose[icos0 - 1]) imin = icos0;
            else imax = icos0;
            icos0 = (imin + imax) / 2;
        }
        if (icos0 == 0) icos0 = 1;
        cos_tmp1 = (cos0 - cose[icos0 - 1]) / (cose[icos0] - cose[icos0 - 1]);
        if (cos_tmp1 < 0.) cos_tmp1 = 0.;
        if (cos_tmp1 > 1.) cos_tmp1 = 1.;
        cos_tmp0 = 1. - cos_tmp1;
        cos_tmp[0]=cos_tmp0;
        cos_tmp[1]=cos_tmp1;
/* 
    // given density, find the corresponding index in logden[]:
        imin = 0;
        imax = ngam;
        igam0 = ngam / 2;
        while ((imax - imin) > 1) {
            if (gamma0 >= gam[igam0 - 1]) imin = igam0;
            else imax = igam0;
            igam0 = (imin + imax) / 2;
        }
        if (igam0 == 0) igam0 = 1;
        gam_tmp1 = (gamma0 - gam[igam0 - 1]) / (gam[igam0] - gam[igam0 - 1]);
        if (gam_tmp1 < 0.) gam_tmp1 = 0.;
        if (gam_tmp1 > 1.) gam_tmp1 = 1.;
        gam_tmp0 = 1. - gam_tmp1; 
 */
    printf("index of abun0 = %d + %e \n",iabun0,abun_tmp1);
    printf("index of xi0 = %d + %e \n",ixi0,xi_tmp1);
    printf("index of logden0 = %d + %e \n",ilogden0,den_tmp1);
    printf("index of cos0 = %d + %e \n",icos0,cos_tmp1);


        flux0 = (double *) malloc(ne_loc *ngam * sizeof(double));

//obtain the emission[gamma,ie] in indice, logxi0, abun0, cos0, logden0
        for(igam=0;igam<ngam;igam++)
            for (ie = 0; ie < ne_loc - 1; ie++) {
                /* y0000 */y[0][0][0][0]=emission[ie + (ne_loc - 1) * (icos0 - 1) + 
                    (ne_loc - 1) * ncose * (ilogden0 - 1)+ 
                    (ne_loc - 1) * ncose * nlogden * (ixi0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * (iabun0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y0001 */y[0][0][0][1]=emission[ie + (ne_loc - 1) * icos0  + 
                    (ne_loc - 1) * ncose * (ilogden0 - 1)+ 
                    (ne_loc - 1) * ncose * nlogden * (ixi0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * (iabun0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y0010 */y[0][0][1][0]=emission[ie + (ne_loc - 1) * (icos0 - 1) + 
                    (ne_loc - 1) * ncose * ilogden0 + 
                    (ne_loc - 1) * ncose * nlogden * (ixi0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * (iabun0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y0011 */y[0][0][1][1]=emission[ie + (ne_loc - 1) * icos0 + 
                    (ne_loc - 1) * ncose * ilogden0 + 
                    (ne_loc - 1) * ncose * nlogden * (ixi0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * (iabun0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y0100 */y[0][1][0][0]=emission[ie + (ne_loc - 1) * (icos0 - 1) + 
                    (ne_loc - 1) * ncose * (ilogden0 - 1)+ 
                    (ne_loc - 1) * ncose * nlogden * ixi0 + 
                    (ne_loc - 1) * ncose * nlogden * nxi * (iabun0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y0101 */y[0][1][0][1]=emission[ie + (ne_loc - 1) * icos0  + 
                    (ne_loc - 1) * ncose * (ilogden0 - 1)+ 
                    (ne_loc - 1) * ncose * nlogden * ixi0  + 
                    (ne_loc - 1) * ncose * nlogden * nxi * (iabun0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y0110 */y[0][1][1][0]=emission[ie + (ne_loc - 1) * (icos0 - 1) + 
                    (ne_loc - 1) * ncose * ilogden0 + 
                    (ne_loc - 1) * ncose * nlogden * ixi0  + 
                    (ne_loc - 1) * ncose * nlogden * nxi * (iabun0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y0111 */y[0][1][1][1]=emission[ie + (ne_loc - 1) * icos0  + 
                    (ne_loc - 1) * ncose * ilogden0 + 
                    (ne_loc - 1) * ncose * nlogden * ixi0 + 
                    (ne_loc - 1) * ncose * nlogden * nxi * (iabun0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y1000 */y[1][0][0][0]=emission[ie + (ne_loc - 1) * (icos0 - 1) + 
                    (ne_loc - 1) * ncose * (ilogden0 - 1)+ 
                    (ne_loc - 1) * ncose * nlogden * (ixi0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * iabun0 + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y1001 */y[1][0][0][1]=emission[ie + (ne_loc - 1) * icos0  + 
                    (ne_loc - 1) * ncose * (ilogden0 - 1)+ 
                    (ne_loc - 1) * ncose * nlogden * (ixi0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * iabun0  + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y1010 */y[1][0][1][0]=emission[ie + (ne_loc - 1) * (icos0 - 1) + 
                    (ne_loc - 1) * ncose * ilogden0 + 
                    (ne_loc - 1) * ncose * nlogden * (ixi0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * iabun0  + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y1011 */y[1][0][1][1]=emission[ie + (ne_loc - 1) * icos0 + 
                    (ne_loc - 1) * ncose * ilogden0 + 
                    (ne_loc - 1) * ncose * nlogden * (ixi0 - 1) + 
                    (ne_loc - 1) * ncose * nlogden * nxi * iabun0 + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y1100 */y[1][1][0][0]=emission[ie + (ne_loc - 1) * (icos0 - 1) + 
                    (ne_loc - 1) * ncose * (ilogden0 - 1)+ 
                    (ne_loc - 1) * ncose * nlogden * ixi0 + 
                    (ne_loc - 1) * ncose * nlogden * nxi * iabun0 + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y1101 */y[1][1][0][1]=emission[ie + (ne_loc - 1) * icos0  + 
                    (ne_loc - 1) * ncose * (ilogden0 - 1)+ 
                    (ne_loc - 1) * ncose * nlogden * ixi0  + 
                    (ne_loc - 1) * ncose * nlogden * nxi * iabun0 + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y1110 */y[1][1][1][0]=emission[ie + (ne_loc - 1) * (icos0 - 1) + 
                    (ne_loc - 1) * ncose * ilogden0 + 
                    (ne_loc - 1) * ncose * nlogden * ixi0  + 
                    (ne_loc - 1) * ncose * nlogden * nxi * iabun0  + 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                /* y1111 */y[1][1][1][1]=emission[ie + (ne_loc - 1) * icos0  + 
                    (ne_loc - 1) * ncose * ilogden0 + 
                    (ne_loc - 1) * ncose * nlogden * ixi0 + 
                    (ne_loc - 1) * ncose * nlogden * nxi * iabun0+ 
                    (ne_loc - 1) * ncose * nlogden * nxi * nabun * igam] ;
                    //printf("%e  %e  %e  %e  \n",y[1][0][0][0],y[0][1][0][0],y[0][0][1][0],y[0][0][0][1]);
            //*0000: first digit for abundence, 2nd for logxi, 3rd for density, 4rd for cose           
                for(i=0;i<2;i++)
                    for(j=0;j<2;j++)
                        for(k=0;k<2;k++)
                            for(m=0;m<2;m++){
                                flux0[ie + ne_loc * igam] += y[i][j][k][m]*abun_tmp[i]*xi_tmp[j]*den_tmp[k]*cos_tmp[m];
                            }
            }
           
        //write_flux(flux0,ne_loc,"flux0.dat");  


/***********************rebin the spectra to the evenly log spaced energies******************************/
        /* nener = 600;
        elow = energy0[0];
        ehigh = energy0[ne_loc-1];

        flux1 = (double *) malloc(nener * ngam *sizeof(double));
        energy2 = (double *) malloc((nener + 1) * sizeof(double));



        energy2[0] = elow / (1. + pow(ehigh / elow, 1. / (nener - 1.))) * 2.;
        for (ie = 1; ie <= nener; ie++) {
            energy2[ie] = energy2[0] * pow(ehigh / elow, ie / (nener - 1.));
        }

        for(igam = 0; igam < ngam; igam++)
            for (ie = 0; ie < nener; ie++) flux1[ie + nener * igam] = 0.;

        ie = 1;
        while (energy0[ie] <= energy2[0]) ie++;
        je = 1;
        quit = 0;
        while ((ie <= (ne_loc - 1)) && (energy0[ie - 1] < energy2[nener])) {
            egap_tmp = energy0[ie] - energy0[ie - 1];
            while (energy2[je - 1] < energy0[ie]) {
            if (energy2[je - 1] < energy0[ie - 1]) etmp = energy0[ie - 1];
            else etmp = energy2[je - 1];
            if (energy2[je] < energy0[ie]) etmp1 = energy2[je];
            else etmp1 = energy0[ie];
            if (etmp1 > etmp)
                for(igam=0;igam<ngam;igam++)
                    flux1[je - 1 + nener * igam] += flux0[ie - 1 + ne_loc * igam] * (etmp1 - etmp) / egap_tmp;
            if (je < nener) je++;
            else {
                quit = 1;
                break;
            }
        }
            if (quit) break;
            if (energy2[je - 1] > energy0[ie]) je--;
            ie++;
            }
        for(igam = 0; igam < ngam; igam++)
            for (ie = 0; ie < nener; ie++) 
                flux1[ie + nener * igam] /= (energy2[ie + 1] - energy2[ie]);
        */
        read_gamma(0);
       
        for(it=0;it<tbin;it++){
            imin = 0;
            imax = ngam;
            igam0 = ngam / 2;
            while ((imax - imin) > 1) {
                if (gam0[it] >= gam[igam0 - 1]) imin = igam0;
                else imax = igam0;
                igam0 = (imin + imax) / 2;
            }
            if (igam0 == 0) igam0 = 1;
            gam_tmp1 = (gam0[it] - gam[igam0 - 1]) / (gam[igam0] - gam[igam0 - 1]);
            if (gam_tmp1 < 0.) gam_tmp1 = 0.;
            if (gam_tmp1 > 1.) gam_tmp1 = 1.;
            gam_tmp0 = 1. - gam_tmp1; 
            printf("gamma index = %d \n",igam0);
            for(ie=0;ie<ne_loc;ie++){
                flux_fin[it][ie] = gam_tmp1*flux0[ie+ne_loc*igam0]+gam_tmp0*flux0[ie+ne_loc*(igam0-1)];
            }
        }
        //write_data(flux0,3000,"flux0.dat");
        write_fluxfin();
        
       return 0;
}






void locate_hdu(fptrp) 
fitsfile *fptrp;
{
    int chdunum=0;
    ffghdn(fptrp,&chdunum);
    printf("The current HDU is HDU  %d \n",chdunum);
}


void write_data(flux,lens,name)
double flux[];
int lens;
char name;
{
    int idx;
    FILE *fp;
    fp=fopen(name,"w");
    for(idx=0;idx<lens;idx++){
        fprintf(fp,"%e  ",flux[idx]);
    }
    fclose(fp);

}/**/


void read_gamma(print)
int print;
{
    FILE *fp_gamma;
   

    fp_gamma = fopen("gamma_ring10_128.dat","r");
    if(fp_gamma == NULL){
            printf("ERROR!\n");
            return;
        }
    rewind(fp_gamma);
    for(i=0;i<tbin;i++){
        fscanf(fp_gamma,"%lf", &gam0[i]);
        gam0[i] *=2;
    }
    printf("gamma0 = %e",gam0[0]);
    if(print) for(i=0;i<tbin;i++)  printf("%e \n",gam0[i]);
    
}

 void write_fluxfin()
{
    FILE *fpfin;
    fpfin=fopen("flux_fin.dat","w");
    for(i=0;i<tbin;i++){
        for(j=0;j<3000;j++){
           fprintf(fpfin,"%e  ",flux_fin[i][j]);
        }
        fprintf(fpfin,"\n");
    }
    fclose(fpfin); 
} 