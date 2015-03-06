//!########################################################################
//!                                                                       #
//! Copyright (C) University of Houston .                     #
//! This file is a finite difference modeling code for 2d acoustic media  #
//!  Usage: Generate multi-shots  reflection   data                                                               #
//! 
//!                                                     		  # 
//!                                                                       #
//!########################################################################
#define _USE_MATH_DEFINES
#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define  MASTER		0

#include "init.h"
#include "output.h"
#include "FM_acoustic.h"

//#include "Reflectiondata.h"

int main( int argc, char *argv[] )
{


  int i,j,m;


  /* Define grids dimensions */
  int dim_w[2];  /* full data dimensions*/
  int nz    = 401;           /* depth sampling nodes   */
  int nx    = 801;	     /* lateral sampling nodes */
  int ne    = 50;            /* absorbing nodes        */
  
     // total grid points number




  // physical grid size
  double dz = 5.0;     // modeling grid spacial step in z
  double dx = 5.0;    // modeling grid spacial step in x

  // time discretization parameter
  int   nt = 1200;        // number of steps
  double dt = 0.001;      // time step

 

 
  /* file for  recording progress */
  FILE* file1;
  char fname1[100];
  int ret;
  sprintf(fname1,"test_2D_acoustic.out");
  file1 = fopen(fname1,"w");
  fprintf(file1, "--!\t                                     \n");
  fprintf(file1, "--!\t2D acoustic modeling research code\n");
  fprintf(file1, "--!\t                                     \n");

  ret = fflush(file1);
   const char * NameWeekDay[] = {"Sunday", "Monday", "Tuesday", "Wenesday", "Thursday", "Friday", "Saturday"};


  time_t timestamp;
  struct tm * tc;
  //struct tm * tf;
    
  timestamp = time(NULL);
  tc = localtime(&timestamp);
  fprintf(file1,"--!\t                                     \n");
  fprintf(file1,"--!\tStart Time : %s, ", NameWeekDay[tc->tm_wday]);
  fprintf(file1,"%02u/%02u/%04u, ", tc->tm_mon+1, tc->tm_mday, 1900 + tc->tm_year);
  fprintf(file1,"Clock : %02uh %02umn %02usec.\n", tc->tm_hour, tc->tm_min, tc->tm_sec);

  ret = fflush(file1);

   // Read main parameter file 
  char file_v0[100],file_sou[100];
  int veloread, nbsou,nbrec;
  double f0;
  params_read(argv[1], &nx, &nz,&dx,&dz,&ne,&nt,&dt,&f0,&veloread, file_v0,&nbsou,&nbrec,file_sou);
  //fprintf(stdout, "%s",file_sou);
 
  int ndims = nx*nz;    
  float* vel  = malloc( ndims*sizeof(double) );  
  if(veloread==0){
  	init_acou_layer( file1,dx,dz, nx,nz,  vel);
  }else{
  	init_acous_vel_read(file1,file_v0,  dx, dz, nx,  nz, vel);
  }
  
  int nz2=2*nz-1;
  int nxe=nx+2*ne;
  int nze=nz2+2*ne;
  int ndims2 = nxe*nze;   
  float* velext  = malloc( ndims2*sizeof(double) );   
  init_acou_vel_ext(file1, dx, dz,  nx,  nz,ne, vel, velext);  


  double* source = malloc(nt*sizeof(double) );
  init_source_ricker_fwps( file1, nt, dt, f0,source); 


   int* psindex = malloc(2*nbsou*sizeof(int));
   double* pscor = malloc(2*nbsou*sizeof(double));
   init_source_cord(file1,  nbsou, nx, nz, ne, file_sou, psindex, pscor);
   
   int nbrec2=(nbrec-1)/2;
   int* prec = malloc(2*nbrec*sizeof(int));
   double* precor = malloc(2*nbrec*sizeof(double));

   int psin[2];

   int ndimt = nbrec*nt; 
   double* P0trace = malloc(ndimt*sizeof(double)); 
   float* P_shotsum = malloc(nbsou*ndimt*sizeof(float));  

   double x,t;


      FILE* file3;  
      char fname3[100];
      FILE* file4;  
      char fname4[100];



    int sousize, soui, souf, nbsouid;
    int numtasks, taskid;

     timestamp = time(NULL);
     tc = localtime(&timestamp);
     fprintf(file1,"--!\t                                     \n");
     fprintf(file1,"--!\tMPI start:");
     fprintf(file1,"%02u/%02u/%04u, ", tc->tm_mon+1, tc->tm_mday, 1900 + tc->tm_year);
     fprintf(file1,"Clock : %02uh %02umn %02usec.\n", tc->tm_hour, tc->tm_min, tc->tm_sec);

     ret = fflush(file1);


	MPI_Init(&argc, &argv);

        MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);  
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);     

	sprintf(fname3,"proc_%d.dat", taskid);
	file3 = fopen(fname3,"w");

     sousize = (nbsou/numtasks)+1; 
     soui = taskid*sousize;  
     souf = (taskid+1)*sousize;     
     if(soui>nbsou) soui = nbsou;
     if(souf>nbsou) souf = nbsou;
     nbsouid = souf-soui;
 
	fprintf (file3, "MPI task %d has started... soui %d souf %d \n", taskid, soui, souf);
        ret = fflush(file3);
  
   for(i=soui;i<souf;i++){ 
		
		psin[0]=psindex[2*i];
		psin[1]=psindex[2*i+1];

		
		for(m=0;m<nbrec;m++){
    			prec[2*m] =  psin[0];
    			prec[2*m+1] =  m-nbrec2+psin[1];
    			precor[2*m] = prec[2*m]*dz;
    			precor[2*m+1] = (prec[2*m+1]-ne)*dx;
   
  		}
               
		
		FM_Acou2D( file3,nx,nz, dx,dz, dt, nt, ne, velext, source, psin,nbrec, prec,P0trace);
		for(j=0;j<ndimt;j++){
			P_shotsum[i*ndimt+j]=(float)P0trace[j];
		}	
    }

   
     int offset, lsize;
     int tag1=2, tag2=1, tag3=3;
     int dest;

     offset = soui*ndimt;
     lsize =(souf-soui)*ndimt;

      
    
     if (taskid > MASTER) {
	 
	  dest = MASTER;
      	  ret=fflush(file3);
	  MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
	  MPI_Send(&P_shotsum[offset], lsize, MPI_DOUBLE, MASTER, tag2, MPI_COMM_WORLD);
     }

     fprintf(file3," after sending \n");
     ret = fflush(file3);

  	
   
     if (taskid == MASTER){

	  
	  for (i=1; i<numtasks; i++) {

	     soui = i*sousize;
	     souf = (i+1)*sousize;     
	     if(soui>nbsou) soui = nbsou;
	     if(souf>nbsou) souf = nbsou;
    		 offset = soui*ndimt;
     		lsize =(souf-soui)*ndimt;


        fprintf(file3," master receive slave %d \n", i);
        ret = fflush(file3);
	MPI_Recv(&offset, 1, MPI_INT, i, tag1, MPI_COMM_WORLD, &status);
        ret = fflush(file3);
        MPI_Recv(&P_shotsum[offset], lsize, MPI_DOUBLE, i, tag2,MPI_COMM_WORLD, &status);
	

         }
     }


    
     if(taskid==MASTER){


   	sprintf(fname4,"reflection.dat");
	file4 = fopen(fname4,"w");
   	   FILE* file6;  
        char fname6[100];
	sprintf(fname6,"reflection");
	file6 = fopen(fname6,"a");
	for (m=0;m<nbsou;m++){
		for (j=0;j<nbrec;j++){
			for (i=0;i<nt;i++){
				t=i*dt;		
				x=(m*nbrec+j)*dx;
				fprintf(file4," %lf  %lf %lf \n", x,t, P_shotsum[m*ndimt+j*nt+i]);
			}
			fprintf(file4," \n");
		}
		fprintf(file4," \n");
	}

	  if( fwrite( P_shotsum, sizeof(float), nbsou*ndimt, file6 ) != (size_t) nbsou*ndimt)
     {
	printf(" Can't output data to file \n");
	exit(1);
     }
      fclose(file6);
    
       
     }

     MPI_Finalize();
    
   
       timestamp = time(NULL);
       tc = localtime(&timestamp);
	  fprintf(file1,"--!\t                                     \n");
	  fprintf(file1,"--!\tModeling finished:");
	  fprintf(file1,"%02u/%02u/%04u, ", tc->tm_mon+1, tc->tm_mday, 1900 + tc->tm_year);
	  fprintf(file1,"Clock : %02uh %02umn %02usec.\n", tc->tm_hour, tc->tm_min, tc->tm_sec);

	  ret = fflush(file1);
	 return (0);


}
