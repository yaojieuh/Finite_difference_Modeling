#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <time.h>

#include "init.h"
#include "FFT.h"
void init_acou_homog(FILE *filefprintf, int *dim_w, int ne, double *vel){
 
 int i;
 //long index;
 int nze = 2*ne + dim_w[0]; 
 int nxe = 2*ne + dim_w[1];
 int ndims = nze*nxe;

 fprintf(filefprintf, "--!\tAcoustic media initialization homgeneous ...\n");
 fprintf(filefprintf, "--!\t                                \n");

  for(i=0; i<ndims; i++) vel[i] = 1500.0;

}

void init_acou_layer(FILE *filefprintf, double dx,double dz, int nx, int nz, float *vel){
 
 int i, j;
 long index;
 double x,z;

// int nzl1= 237, nzl2=504;
 int nzl1= 30, nzl2=55;
 int ndims=nx*nz;
 

 fprintf(filefprintf, "--!\tAcoustic media initialization three layers ...\n");
 fprintf(filefprintf, "--!\t                                \n");
	FILE* file2;  
  	char fname2[100];

  	sprintf(fname2,"Velocitymodel");
  	file2= fopen( fname2, "w" );
	  for(j=0; j<nx; j++){
		x=j*dx;
	      	for(i=0; i<nz; i++) {
			z=i*dz;
			index= j*nz+i;
			if (i<nzl1){
				vel[index]=1500;
			}
		 	else if (i<nzl2){
                  		vel[index]=2000;
                 
	      		}
	     		else{
                  		vel[index]=1500;
                  	 }
	 		
	 	 }		
	}
          if( fwrite( vel, sizeof(float), ndims, file2 ) != (size_t) ndims)
     {
	printf(" Can't output data to file \n");
	exit(1);
     }    
	fclose(file2);
}

void init_acous_vel_read(FILE *filefprintf,char *file_v0, double dx,double dz, int nx, int nz, float *vel){

 int ndims=nx*nz;
 

 FILE* file4;
 char fname4[100];
 sprintf(fname4,file_v0);
 file4 = fopen(fname4, "rb");
 size_t size_e;
  if ( !file4 )
  {
    fprintf(stderr, "\n Problem opening file %s ", fname4);
    exit(EXIT_FAILURE);
  }
  size_e = fread( vel, sizeof(float), ndims, file4 );
  if( size_e != (size_t) ndims )
  {
    fprintf(stderr, "\n Problem reading file %s ", fname4);
    exit(EXIT_FAILURE);
  }
  fclose(file4);

}


void init_acou_vel_ext(FILE *filefprintf, double dx,double dz, int nx, int nz,int ne, float *vel, float *velext){
 
 int ix, iz;
 int index, index2;
 double x,z;
 int nz2=2*nz-1;
 int nxe=nx+2*ne;
 int nze=nz2+2*ne;
 int ndims=nx*nz;
 int ndims2=nxe*nze;
 

 fprintf(filefprintf, "--!\tAcoustic media extention ...\n");
 fprintf(filefprintf, "--!\t                                \n");

 FILE* file2;  
 char fname2[100];

 sprintf(fname2,"Velocitymodel_ext");
 file2= fopen( fname2, "w" );
	  for(ix=0; ix<nxe; ix++){
	      	for(iz=0; iz<nze; iz++) {
			index= ix*nze+iz;
			if (iz<ne+nz-1){
				velext[index]=vel[0];
			}
		 	else if (iz<ne+2*nz-1){
                  		if (ix<ne){
					velext[index]=vel[iz-(ne+nz-1)];
				}
				else if(ix<ne+nx){
					index2=(ix-ne)*nz+iz-(ne+nz-1);
					velext[index]=vel[index2];
				}
				else{
					index2=(nx-1)*nz+iz-(ne+nz-1);
					velext[index]=vel[index2];
				}
                 
	      		}
	     		else{
                  		if (ix<ne){
					velext[index]=vel[nz-1];
				}
				else if(ix<ne+nx){
					index2=(ix-ne)*nz+nz-1;
					velext[index]=vel[index2];
				}
				else{
					index2=(nx-1)*nz+nz-1;
					velext[index]=vel[index2];
				}
                  	 }
	 		
	 	 }		
	}
          if( fwrite( velext, sizeof(float), ndims2, file2 ) != (size_t) ndims2)
     {
	printf(" Can't output data to file \n");
	exit(1);
     }    
	fclose(file2);
}

void init_source_ricker_fwps(FILE *filefprintf, int nt, double dt, double f0, double *source){

  int i=0; 
  double at=0.0, arg=0.0;
  double sig=1.5, gam=8.0, tau=1.0;
  double pi=3.14159265358979323846;
  double fmax=f0;
  FILE* file2;  
  char fname2[100];
  fprintf(filefprintf, "--!\tSource wavelet initialization  ...\n");
  fprintf(filefprintf, "--!\t                                   \n");

  sprintf(fname2,"Ricker_fwps_%.2lf_Hz.dat",fmax);
  file2= fopen( fname2, "w" );

  sig *= fmax;
  for(i=0; i<nt; i++ ) {

         at = i * dt;
         arg = -1.0*sqrt(2.0/pi)*sig*gam;
         arg = arg*(sig-2.0*sig*gam*(sig*at-tau)*(sig*at-tau));
         source[i] = 0.05*arg*exp(-gam*(sig*at-tau)*(sig*at-tau));
         fprintf(file2, " %lf %lf \n", i*dt, source[i]);

     }
		              
     	fclose(file2);
	

}
void init_source_gauss(FILE *filefprintf, int nt, double dt, double dw, double *source, double *sw, int nfmax){

  int i=0; 
  double at=0.0,arg;
  double sig=0;
  double pi=3.14159265358979323846;
  double fmax=20;
  FILE* file2;  
  char fname2[100];
  int nf;
  fprintf(filefprintf, "--!\tSource wavelet initialization  ...\n");
  fprintf(filefprintf, "--!\t                                   \n");

  sprintf(fname2,"Ricker_fwps_%.2lf_Hz.dat",fmax);
  file2= fopen( fname2, "w" );

  sig = 100*pi*pi;
  for(i=0; i<nt; i++ ) {
         at = (i * dt-0.2);
         arg = -at*at*sig;
         source[i] = 20*at*exp(arg)*sig*sqrt(pi);
         fprintf(file2, " %lf %lf \n", i*dt, source[i]);

     }
		              
     	fclose(file2);
	FILE* file3;  
     	char fname3[100];
	sprintf(fname3,"Rickerfren_fwps_%.2lf_Hz.dat",fmax);
	file3 = fopen(fname3,"w");
	double t, aw=32.0, ac,agr;
	int iw;
    	 for(iw=0;iw<nfmax;iw++){

            aw = iw*dw + 0.001;
	    sw[iw*3]=aw;
            at = 0.0;
            arg = aw*at;

            ac = 0.5; 
	    sw[3*iw+1]  =      ac*dt*cos(arg)*source[0];
	    sw[3*iw+2] = -1.0*ac*dt*sin(arg)*source[0];

            ac = 1.0;
	    for(i=1;i<nt-1;i++){
                at = i*dt;
                arg = aw*at;
		sw[3*iw+1] += ac*dt*cos(arg)*source[i];
		sw[3*iw+2] -= ac*dt*sin(arg)*source[i];
	    }

                at = (nt-1)*dt;
                ac = 0.5;
                arg = aw*at;
		sw[3*iw+1]  += ac*dt*cos(arg)*source[nt-1];
		sw[3*iw+2] -= ac*dt*sin(arg)*source[nt-1];
		//sw[3*iw+1]  = 0;
		//sw[3*iw+2] *= 2;
    
         fprintf(file3," %lf %lf %lf  %lf\n", aw, sw[3*iw+1], sw[3*iw+2], -aw*exp(-aw*aw/4000));
	 }
	fclose(file3);

}

void init_source_cord(FILE *file1, int nbsou,int nx,int nz,int ne, char *file_sou, int *psindex, double * pscor){

 
   FILE* file2;
   char fname2[100];
   sprintf(fname2,file_sou);
   file2 = fopen(fname2,"r");
   fprintf(file1," Read source coordinate from %s \n",fname2);	
   int i,readtemp;
   int* readdat = malloc((2*nbsou+1)*sizeof(int));
	while (!feof(file2))
	{
		fscanf(file2, "%d",&readtemp);
		readdat[i]=readtemp;
                
		i++;
	}
	for (i=0;i<nbsou;i++){
		psindex[2*i]=readdat[2*i]+(nz-1)+ne;
		psindex[2*i+1]=readdat[2*i+1]+ne;
		
	}
      fclose(file2);
 



}
