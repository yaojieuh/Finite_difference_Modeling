#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <time.h>

//#include "output.h"
#include  "operation_complex.h"
void comput_rec(int it, int nbrec, int *prec, double *u3, int nx,int nz, int ne, double dx,double dz, double *recs){

   int i;
   int ix, iz; 
   int i1;
   
   int nze =2*nz-1 + 2*ne;
   int nxe = nx+ 2*ne;

   for(i=0;i<nbrec;i++){

      iz = prec[i*2]; 
      ix = prec[i*2+1]; 

      i1 = iz  + (ix)*nze;

      recs[i+it*nbrec] = u3[i1]/dx/dz;  //p-wave value at receiver, position,   for time step =it 

    }  

}


// Output 2D p-wave snapshots 
void Output_snapshots_xz( int ne, int it, double x0, int *dim_w, double *dh, double *u2, double *uv){

     int i, j, k, nxe, nze;
     long index, index1, ndims; 
     nxe = dim_w[1]+2*ne; 		
     nze = dim_w[0]+2*ne; 
     ndims = nxe*nze;

     FILE* file1;  
     char fname1[100];
		sprintf(fname1,"snapshot_source_%.2lf_xz_%d.dat",x0,it);
		file1 = fopen(fname1,"w");

      for(i=0;i<dim_w[0];i++){
         for(j=0;j<dim_w[1];j++){
	  index = i + ne + (j+ne)*nze;
	  index1 = i + j*dim_w[0];
	  uv[index1] = u2[index];
          fprintf(file1," %lf %lf %lf \n",dh[1]*j,-dh[0]*i,uv[index1]);
           }
           fprintf(file1,"  \n");
      }
 
}




void Output_rec_trac_acou( int nbrec, int *prec, int ne, double x0, double dx, double dz,int it, double dt,  double *recs){

     int i=0, j=0;
     double x; 
     		FILE* file1;  
     		char fname1[100];
		sprintf(fname1,"map_rec_%.2lf.dat",x0);
		file1 = fopen(fname1,"w");
 		

         for(i=0;i<nbrec;i++){
		x=(prec[2*i+1])*dx;
		for(j=0;j<it;j++){
		  fprintf(file1," %f %f %f\n", x,(j+1)*dt,  recs[j*nbrec+i]);
		}
              fprintf(file1,"\n");
     	 }

}

