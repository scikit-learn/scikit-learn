
#if defined(__LP64__) /* In LP64 match sizes with the 32 bit ABI */
typedef int 		ml_int;
#else
typedef long int 	ml_int;
#endif



/* lars.c */
void lars_fit(int nfeatures, int nsamples, double *X, double *b, double *beta,double *lambdas,  int *row, int *col, double *L, int stop);

/* from cholesky.c */
int cholesky_update(double *R, int ldr, int p, double *X, double *Z, int ldz,
                    int nz, double *Y, double *rho, double *C, double *S);

void givens_rot(double *da, double *db, double *dc, double *ds);

double cholesky_logdet (double *L, int n);


