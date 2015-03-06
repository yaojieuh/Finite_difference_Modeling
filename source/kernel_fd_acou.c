#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "kernel_fd_acou.h"

void Kernel_fd_acou( int nx, int nz, int ne, double dx, double dz, double dt, double *ax, double *bx, double *az, double *bz, float *v, double *u3, double *u1, double *u2,double *psix1, double *psix2, double *psiz1, double *psiz2){

  long i, index, j, k;
  int nze =2*nz-1 + 2*ne;
  int nxe = nx+ 2*ne;
  long ndims = nxe*nze;
  double pzm1,pzm2,pzm3,pzm4;
  double pzp1,pzp2,pzp3,pzp4;
  double pxm1,pxm2,pxm3,pxm4;
  double pxp1,pxp2,pxp3,pxp4;
  double hx, hz;
  int dh[2];
  dh[0]=dz, dh[1]=dx;

  /* stencil coeffs */
  double a0 = -2.84722222, a1 = 1.6, a2 = -0.2;
  double a3 = 0.0253968254, a4 = -0.00178571429;
         
  /* intervals */
  double dtpdz = dt/dh[0];
  double dtpdx = dt/dh[1];
	
  double odz2 = dtpdz*dtpdz;
  double odx2 = dtpdx*dtpdx;

  /* coeffs, b, 8th order */
  /* z coeffs */
  double b0z  = a0*odz2;
  double b1z = a1 * odz2;
  double b2z = a2 * odz2;
  double b3z = a3 * odz2;
  double b4z = a4 * odz2;

  /* coeffs, c, 2nd order */
  /*  z coeffs */
  double c0z  = -2 * odz2;
  double c1z =  1 * odz2;

  /* coeffs, b, 8th order */
  /* x coeffs */
  double b0x  = a0*odx2;
  double b1x = a1 * odx2;
  double b2x = a2 * odx2;
  double b3x = a3 * odx2;
  double b4x = a4 * odx2;

  /* coeffs, c, 2nd order */
  /*  x coeffs */
  double c0x  = -2 * odx2;
  double c1x =  1 * odx2;

//////////////////////////////////////////////////////////////////
// compute d/dx2+d/dz2 pn
  for(j=1; j<ne; j++) {

	     for(i=1; i<4; i++) {

		/* stencil z */
		index = i-1 + j*nze;
		pzm1 = u2[index];
		index = i+1 + j*nze;
		pzp1 = u2[index];

		/* stencil x */
		index = i + (j-1)*nze;
		pxm1 = u2[index];
		index = i + (j+1)*nze;
		pxp1 = u2[index];

		index = i + j*nze;

		hx = c0x*u2[index] + c1x*(pxm1 + pxp1);
		hz = c0z*u2[index] + c1z*(pzm1 + pzp1);

		u3[index] = hx + hz; 
		psix2[index]=bx[j]*psix1[index]+ax[j]*(u2[ i + (j+1)*nze]-u2[ i + (j)*nze])/dh[1];
		psiz2[index]=bz[i]*psiz1[index]+ax[i]*(u2[i+1 + j*nze]-u2[i + (j)*nze])/dh[0];
		u3[index]+=((psix2[index]-psix2[i + (j-1)*nze])/dh[1]+(psiz2[index]-psiz2[i-1 + j*nze])/dh[0])*dt*dt;

	     }

	     for(i=4; i<nze-4; i++) {

		/* stencil z */
		index = i-4 + j*nze;
		pzm4 = u2[index];
		index++;
		pzm3 = u2[index];
		index++;
		pzm2 = u2[index];
		index++;
		pzm1 = u2[index];

		index += 2;
		pzp1 = u2[index];
		index++;
		pzp2 = u2[index];
		index++;
		pzp3 = u2[index];
		index++;
		pzp4 = u2[index];
	  
		/* stencil x */
		index = i + (j-1)*nze ;
		pxm1 = u2[index];
		index = i + (j+1)*nze ;
		pxp1 = u2[index];

		index = i + j*nze ;

		hx = c0x*u2[index] + c1x*(pxm1 + pxp1);

		hz = b4z*(pzm4+pzp4);
		hz += b3z*(pzm3+pzp3);
		hz += b2z*(pzm2+pzp2);
		hz += b1z*(pzm1+pzp1);
		hz += b0z*u2[index];

		u3[index] = hx + hz; 

		psix2[index]=bx[j]*psix1[index]+ax[j]*(u2[ i + (j+1)*nze]-u2[ i + (j)*nze])/dh[1];
		psiz2[index]=bz[i]*psiz1[index]+ax[i]*(u2[i+1 + j*nze]-u2[i + (j)*nze])/dh[0];
		u3[index]+=((psix2[index]-psix2[i + (j-1)*nze])/dh[1]+(psiz2[index]-psiz2[i-1 + j*nze])/dh[0])*dt*dt;

	     }

	     for(i=nze-4; i<nze-1; i++) {


		/* stencil z */
		index = i-1 + j*nze;
		pzm1 = u2[index];
		index = i+1 + j*nze;
		pzp1 = u2[index];

		/* stencil x */
		index = i + (j-1)*nze;
		pxm1 = u2[index];
		index = i + (j+1)*nze;
		pxp1 = u2[index];

		index = i + j*nze;

		hx = c0x*u2[index] + c1x*(pxm1 + pxp1);
		hz = c0z*u2[index] + c1z*(pzm1 + pzp1);

		u3[index] = hx + hz; 
		psix2[index]=bx[j]*psix1[index]+ax[j]*(u2[ i + (j+1)*nze]-u2[ i + (j)*nze])/dh[1];
		psiz2[index]=bz[i]*psiz1[index]+ax[i]*(u2[i+1 + j*nze]-u2[i + (j)*nze])/dh[0];
		u3[index]+=((psix2[index]-psix2[i + (j-1)*nze])/dh[1]+(psiz2[index]-psiz2[i-1 + j*nze])/dh[0])*dt*dt;

	     }

	  }

	
	  
	  for(j=ne; j<nxe-ne; j++) {

	     for(i=1; i<ne; i++) {
		/* stencil z */
		index = i-1 + j*nze;
		pzm1 = u2[index];
		index = i+1 + j*nze;
		pzp1 = u2[index];

		/* stencil x */
		index = i + (j-1)*nze;
		pxm1 = u2[index];
		index = i + (j+1)*nze;
		pxp1 = u2[index];

		index = i + j*nze;

		hx = c0x*u2[index] + c1x*(pxm1 + pxp1);
		hz = c0z*u2[index] + c1z*(pzm1 + pzp1);

		u3[index] = hx + hz; 
		
		psix2[index]=bx[j]*psix1[index]+ax[j]*(u2[ i + (j+1)*nze]-u2[ i + (j)*nze])/dh[1];
		psiz2[index]=bz[i]*psiz1[index]+ax[i]*(u2[i+1 + j*nze]-u2[i + (j)*nze])/dh[0];
		u3[index]+=((psix2[index]-psix2[i + (j-1)*nze])/dh[1]+(psiz2[index]-psiz2[i-1 + (j)*nze])/dh[0])*dt*dt;
	     }


             // CENTER
	     for(i=ne; i<nze-ne; i++) {

		/* stencil z */
		index = i-4 + j*nze;
		pzm4 = u2[index];
		index++;
		pzm3 = u2[index];
		index++;
		pzm2 = u2[index];
		index++;
		pzm1 = u2[index];

		index += 2;
		pzp1 = u2[index];
		index++;
		pzp2 = u2[index];
		index++;
		pzp3 = u2[index];
		index++;
		pzp4 = u2[index];
	  
		/* stencil x */
		index = i + (j-4)*nze;
		pxm4 = u2[index];
		index += nze;
		pxm3 = u2[index];
		index += nze;
		pxm2 = u2[index];
		index += nze;
		pxm1 = u2[index];

		index = i + (j+1)*nze;
		pxp1 = u2[index];
		index += nze;
		pxp2 = u2[index];
		index += nze;
		pxp3 = u2[index];
		index += nze;
		pxp4 = u2[index];

		index = i + j*nze;
                
                // 8th order stencil 
		hx = b4x*(pxm4+pxp4);
		hx += b3x*(pxm3+pxp3);
		hx += b2x*(pxm2+pxp2);
		hx += b1x*(pxm1+pxp1);
		hx += b0x*u2[index];

		hz = b4z*(pzm4+pzp4);
		hz += b3z*(pzm3+pzp3);
		hz += b2z*(pzm2+pzp2);
		hz += b1z*(pzm1+pzp1);
		hz += b0z*u2[index];

		u3[index] = hx + hz; 

	     }

	     for(i=nze-ne; i<nze-1; i++) {
	/* stencil z */
		index = i-1 + j*nze;
		pzm1 = u2[index];
		index = i+1 + j*nze;
		pzp1 = u2[index];

		/* stencil x */
		index = i + (j-1)*nze;
		pxm1 = u2[index];
		index = i + (j+1)*nze;
		pxp1 = u2[index];

		index = i + j*nze;

		hx = c0x*u2[index] + c1x*(pxm1 + pxp1);
		hz = c0z*u2[index] + c1z*(pzm1 + pzp1);;

		u3[index] = hx + hz; 
		psix2[index]=bx[j]*psix1[index]+ax[j]*(u2[ i + (j+1)*nze]-u2[ i + (j)*nze])/dh[1];
		psiz2[index]=bz[i]*psiz1[index]+ax[i]*(u2[i+1 + j*nze]-u2[i + (j)*nze])/dh[0];
		u3[index]+=((psix2[index]-psix2[i + (j-1)*nze])/dh[1]+(psiz2[index]-psiz2[i-1 + (j)*nze])/dh[0])*dt*dt;

	     }

	  }

	  for(j=nxe-ne; j<nxe-1; j++) {


	     for(i=1; i<nze-1; i++) {


		/* stencil z */
		index = i-1 + j*nze;
		pzm1 = u2[index];
		index = i+1 + j*nze;
		pzp1 = u2[index];

		/* stencil x */
		index = i + (j-1)*nze;
		pxm1 = u2[index];
		index = i + (j+1)*nze;
		pxp1 = u2[index];

		index = i + j*nze;

		hx = c0x*u2[index] + c1x*(pxm1 + pxp1);
		hz = c0z*u2[index] + c1z*(pzm1 + pzp1);

		u3[index] = hx + hz; 

		psix2[index]=bx[j]*psix1[index]+ax[j]*(u2[ i + (j+1)*nze]-u2[ i + (j)*nze])/dh[1];
		psiz2[index]=bz[i]*psiz1[index]+ax[i]*(u2[i+1 + j*nze]-u2[i + (j)*nze])/dh[0];
		u3[index]+=((psix2[index]-psix2[i + (j-1)*nze])/dh[1]+(psiz2[index]-psiz2[i-1 + (j)*nze])/dh[0])*dt*dt;

	     }

	  }

//    pn+1 = 2 pn - pn-1 +c^2*dt^2*lap pn
  	for(i=0; i<ndims; i++) u3[i] = 2*u2[i] - u1[i] + v[i]*v[i]*u3[i];

}

void Intro_source_acou_fwps( int it,int nx, int nz, int ne,  int *ps, double dt, double *source, float *v, double *u3){
 
   double rr=0.0;
   //int i, j=0, k, ix, iz, iy;
   int ix, iz; 
   int index;
  // double pi=3.14159265358979323846;

   iz = ps[0] ;
   ix = ps[1] ;

   int nze =2*nz-1 + 2*ne;
   int nxe = nx+ 2*ne;

 

 //	rr = sqrt((j-ix)*(j-ix)+(i-iz)*(i-iz)); 
 	index = iz + ix*nze; 
	rr  = dt*dt*v[index]*v[index]*source[it];  
 	u3[index] = u3[index] - rr; 

}

