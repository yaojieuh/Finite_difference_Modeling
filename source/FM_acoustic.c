#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <time.h>

#include "FM_acoustic.h"
#include "kernel_fd_acou.h"
#include "output.h"
//#include "Pscatter.h"

void FM_Acou2D( FILE *file1, int nx, int nz, double dx, double dz, double dt, int nt, int ne, float *v, double *source, int *psou, int nbrec, int *prec, double *P0trace){
  double x0=(psou[1]-ne)*dx;
  int nze =2*nz-1 + 2*ne;
  int nxe = nx+ 2*ne;
  int ndims = nxe*nze;
  long i=0,j=0;
  int it=0; 
  double t=0.0;
  double* recs = malloc(nt*nbrec*sizeof(double));
  double pi = 3.14159265358979323846; 
  double cst = 1/sqrt(2.0*pi);
  // p wavefields
  double* u1 = malloc( ndims*sizeof(double) );
  double* u2 = malloc( ndims*sizeof(double) );
  double* u3 = malloc( ndims*sizeof(double) );
  double kappax=1, kappaz=1, dampmax=sqrt(200);
  double* psix1 = malloc( ndims*sizeof(double) );
  double* psix2 = malloc( ndims*sizeof(double) );
  double* psiz1 = malloc( ndims*sizeof(double) );
  double* psiz2 = malloc( ndims*sizeof(double) );
  double* dampx = malloc( nxe*sizeof(double) );
  double* alphax = malloc( nxe*sizeof(double) );
  double* ax = malloc( nxe*sizeof(double) );
  double* bx = malloc( nxe*sizeof(double) );
  double* dampz = malloc( nze*sizeof(double) );
  double* alphaz = malloc( nze*sizeof(double) );
  double* az = malloc( nze*sizeof(double) );
  double* bz = malloc( nze*sizeof(double) );

 


   int ret;
  for (i=0;i<nxe;i++){
	if(i<=ne){
		dampx[i]=(dampmax*(i-ne)/ne)*(dampmax*(i-ne)/ne);
		alphax[i]=0.01*i;
	}
	else if (i>nx+ne){
		dampx[i]=(dampmax*(i-nx-ne)/ne)*(dampmax*(i-nx-ne)/ne);
		alphax[i]=0.01*(nxe-i);
	}
	else {
		dampx[i]=0;
		alphax[i]=0.0001;
	}
	bx[i]=exp(-(dampx[i]/kappax+alphax[i])*dt);
	ax[i]=dampx[i]*(bx[i]-1)/kappax/(dampx[i]+kappax*alphax[i]);	
       // fprintf( file1,"--!%lf_%lf  \n",ax[i],bx[i]);
  }
  for (i=0;i<nze;i++){
	if(i<=ne){
		dampz[i]=(dampmax*(i-ne)/ne)*(dampmax*(i-ne)/ne);
		alphaz[i]=0.1*i;
	}
	else if (i>2*nz-1+ne){
		dampz[i]=(dampmax*(i-(2*nz-1)-ne)/ne)*(dampmax*(i-(2*nz-1)-ne)/ne);
		alphaz[i]=0.1*(nze-i);
	}
	else {
		dampz[i]=0;
		alphaz[i]=0.0001;
	}
	bz[i]=exp(-(dampz[i]/kappaz+alphaz[i])*dt);
	az[i]=dampz[i]*(bz[i]-1)/kappaz/(dampz[i]+kappaz*alphaz[i]);
  }		



 


  for(i=0;i<ndims;i++){
	u1[i]  = 0.0;
	u2[i]  = 0.0;
	u3[i]  = 0.0;
	psix1[i] = 0;
	psix2[i] = 0;
	psiz1[i] = 0;
	psiz2[i] = 0;
  }

 
  for(it=0; it<nt;it++) {


    if((it%100)==0) fprintf( file1,"--!\t%lf Propagation Iteration %d  ...\n",x0,it);
    if((it%500)==0) ret = fflush( file1);

    //Kernel_fft_acou( dimw, ne, dh, dt, v, u3, u1, u2);
    Kernel_fd_acou( nx,nz, ne, dx,dz, dt,ax,bx,az,bz, v, u3, u1, u2,psix1,psix2,psiz1,psiz2);


    // Source Wavelet introduction

    Intro_source_acou_fwps( it, nx,nz, ne, psou, dt, source, v, u3);


    // Computing wave amplitude at receivers  
    comput_rec( it, nbrec, prec, u3, nx,nz, ne, dx,dz, recs);
 
    // Update for next iteration
    for(i=0;i<ndims;i++){

	t  = u1[i];
	u1[i] = u2[i];
	u2[i] = u3[i];
	u3[i] = t;
	psix1[i]=psix2[i];
	psiz1[i]=psiz2[i];
	psix2[i]=0;
	psiz2[i]=0;

    }

   

    //if(((it%200)==0)&&(it>0)){
      //   Output 2D p-wave snapshots 
    //   Output_snapshots_xz( ne, it, x0,dimw, dh, u2, uv);
    //}

  
  }


  // Writing wave amplitude file at receivers
 	for(i=0;i<nbrec;i++){
		P0trace[i]=0;
		for(j=0;j<nt;j++){
			P0trace[i*nt+j]=recs[j*nbrec+i];
		}
	}
	if((((psou[1]-ne)%100)==0)){ 
  		Output_rec_trac_acou( nbrec,prec, ne, x0, dx,dz, nt, dt, recs); 	
	}		
	
  free(u1);
  free(u2);
  free(u3);
  free(recs);
  free(psix1); 
  free(psix2);
  free(psiz1);
  free(psiz2);
  free(dampx); 
  free(alphax);
  free(ax);
  free(bx); 
  free(dampz); 
  free(alphaz);
  free(az);
  free(bz); 
  
}

