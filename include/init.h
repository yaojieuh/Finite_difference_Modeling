void init_acou_layer(FILE *filefprintf, double dx,double dz, int nx, int nz, float *vel);
void init_acou_homog(FILE *,  int *, int, double *);
void init_source_ricker_fwps(FILE *filefprintf, int nt, double dt, double f0, double *source);
void init_source_gauss(FILE *filefprintf, int nt, double dt, double dw, double *source, double *sw, int nfmax);
void init_acous_vel_read(FILE *filefprintf,char *file_v0, double dx,double dz, int nx, int nz, float *vel);
void init_acou_vel_ext(FILE *filefprintf, double dx,double dz, int nx, int nz,int ne, float *vel, float *velext);
void init_source_cord(FILE *file1, int nbsou,int nx,int nz,int ne, char *file_sou, int *psindex, double * pscor);
