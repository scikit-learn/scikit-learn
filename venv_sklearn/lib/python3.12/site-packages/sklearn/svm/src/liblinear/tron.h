#ifndef _TRON_H
#define _TRON_H

#include "_cython_blas_helpers.h"

class function
{
public:
	virtual double fun(double *w) = 0 ;
	virtual void grad(double *w, double *g) = 0 ;
	virtual void Hv(double *s, double *Hs) = 0 ;

	virtual int get_nr_variable(void) = 0 ;
	virtual ~function(void){}
};

class TRON
{
public:
	TRON(const function *fun_obj, double eps = 0.1, int max_iter = 1000, BlasFunctions *blas = 0);
	~TRON();

	int tron(double *w);
	void set_print_string(void (*i_print) (const char *buf));

private:
	int trcg(double delta, double *g, double *s, double *r);
	double norm_inf(int n, double *x);

	double eps;
	int max_iter;
	function *fun_obj;
	BlasFunctions *blas;
	void info(const char *fmt,...);
	void (*tron_print_string)(const char *buf);
};
#endif
