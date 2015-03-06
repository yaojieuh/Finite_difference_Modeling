#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "operation_complex.h"


// mutiple two complex number
void prodcomp( double *z1, double *z2, double *zres ){

     // real part
     zres[0] = z1[0]*z2[0]-z1[1]*z2[1]; 
     // imaginary part  
     zres[1] = z1[1]*z2[0]+z1[0]*z2[1]; 

}

// division two complex numbers
void divcomp( double *z1, double *z2, double *zres ){

     double mod2z2 = 0, zresr, zresi; 

     // num z1 z2*
     zresr = z1[0]*z2[0]+z1[1]*z2[1]; 
     zresi = z1[1]*z2[0]-z1[0]*z2[1]; 

     // denom
     mod2z2 = z2[0]*z2[0]+z2[1]*z2[1];

     if(fabs(mod2z2)>0.000000000001){
        zres[0] = zresr / mod2z2; 
        zres[1] = zresi / mod2z2; 
     }else{
        fprintf(stdout," ERROR DIVISION BY ZERO \n");
        exit(0);
     } 

}
