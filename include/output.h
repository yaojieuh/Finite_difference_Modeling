void comput_rec(int it, int nbrec, int *prec, double *u3, int nx,int nz, int ne, double dx,double dz, double *recs);
void Output_snapshots_xz( int ne, int it, double x0, int *dim_w, double *dh, double *u2, double *uv);
void Output_rec_trac_acou( int nbrec, int *prec, int ne, double x0, double dx, double dz,int it, double dt,  double *recs);
void Output_rec_trac_fren( int nbrec, int *prec, int ne,  double x0, double *dh,int nt, double dt, double *recs, double *Pwr, double *Pwi);
