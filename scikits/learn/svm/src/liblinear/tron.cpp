#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <vector>
#include "tron.h"

extern "C" {

extern double dnrm2_(int *, double *, int *);
extern double ddot_(int *, double *, int *, double *, int *);
extern int daxpy_(int *, double *, double *, int *, double *, int *);
extern int dscal_(int *, double *, double *, int *);

}

static void default_print(const char *buf)
{
	fputs(buf,stdout);
	fflush(stdout);
}

void TRON::info(const char *fmt,...)
{
	// char buf[BUFSIZ];
	// va_list ap;
	// va_start(ap,fmt);
	// vsprintf(buf,fmt,ap);
	// va_end(ap);
	// (*tron_print_string)(buf);
}

TRON::TRON(function &fun_obj, double eps, int max_iter)
 : fun_obj(fun_obj),
   eps(eps),
   max_iter(max_iter),
   tron_print_string(&default_print)
{
}

TRON::~TRON()
{
}

void TRON::tron(double *w)
{
	// Parameters for updating the iterates.
	double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;

	// Parameters for updating the trust region size delta.
	double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

	int n = fun_obj.get_nr_variable();
	int i, cg_iter;
	double delta, snorm, one=1.0;
	double alpha, f, fnew, prered, actred, gs;
	int search = 1, iter = 1, inc = 1;

	std::vector<double> v_s(n),
						v_r(n),
						v_w_new(n),
						v_g(n);
	double *s = &v_s[0];
	double *r = &v_r[0];
	double *w_new = &v_w_new[0];
	double *g = &v_g[0];

	for (i=0; i<n; i++)
		w[i] = 0;

	f = fun_obj.fun(w);
	fun_obj.grad(w, g);
	delta = dnrm2_(&n, g, &inc);
	double gnorm1 = delta;
	double gnorm = gnorm1;

	if (gnorm <= eps*gnorm1)
		search = 0;

	iter = 1;

	while (iter <= max_iter && search)
	{
		cg_iter = trcg(delta, g, s, r);

		std::copy(w, w+n, w_new);
		daxpy_(&n, &one, s, &inc, w_new, &inc);

		gs = ddot_(&n, g, &inc, s, &inc);
		prered = -0.5*(gs-ddot_(&n, s, &inc, r, &inc));
		fnew = fun_obj.fun(w_new);

		// Compute the actual reduction.
	        actred = f - fnew;

		// On the first iteration, adjust the initial step bound.
		snorm = dnrm2_(&n, s, &inc);
		if (iter == 1)
			delta = std::min(delta, snorm);

		// Compute prediction alpha*snorm of the step.
		if (fnew - f - gs <= 0)
			alpha = sigma3;
		else
			alpha = std::max(sigma1, -0.5*(gs/(fnew - f - gs)));

		// Update the trust region bound according to the ratio of actual to predicted reduction.
		if (actred < eta0*prered)
			delta = std::min(std::max(alpha, sigma1)*snorm, sigma2*delta);
		else if (actred < eta1*prered)
			delta = std::max(sigma1*delta, std::min(alpha*snorm, sigma2*delta));
		else if (actred < eta2*prered)
			delta = std::max(sigma1*delta, std::min(alpha*snorm, sigma3*delta));
		else
			delta = std::max(delta, std::min(alpha*snorm, sigma3*delta));

		info("iter %2d act %5.3e pre %5.3e delta %5.3e f %5.3e |g| %5.3e CG %3d\n", iter, actred, prered, delta, f, gnorm, cg_iter);

		if (actred > eta0*prered)
		{
			iter++;
			memcpy(w, w_new, sizeof(double)*n);
			f = fnew;
		        fun_obj.grad(w, g);

			gnorm = dnrm2_(&n, g, &inc);
			if (gnorm <= eps*gnorm1)
				break;
		}
		if (f < -1.0e+32)
		{
			info("warning: f < -1.0e+32\n");
			break;
		}
		if (fabs(actred) <= 0 && prered <= 0)
		{
			info("warning: actred and prered <= 0\n");
			break;
		}
		if (fabs(actred) <= 1.0e-12*fabs(f) &&
		    fabs(prered) <= 1.0e-12*fabs(f))
		{
			info("warning: actred and prered too small\n");
			break;
		}
	}
}

int TRON::trcg(double delta, double *g, double *s, double *r)
{
	int i, inc = 1;
	int n = fun_obj.get_nr_variable();
	double one = 1;
	std::vector<double> v_d(n),
						v_Hd(n);
	double *d = &v_d[0];
	double *Hd = &v_Hd[0];
	double rTr, rnewTrnew, alpha, beta, cgtol;

	for (i=0; i<n; i++)
	{
		s[i] = 0;
		r[i] = -g[i];
		d[i] = r[i];
	}
	cgtol = 0.1*dnrm2_(&n, g, &inc);

	int cg_iter = 0;
	rTr = ddot_(&n, r, &inc, r, &inc);
	while (1)
	{
		if (dnrm2_(&n, r, &inc) <= cgtol)
			break;
		cg_iter++;
		fun_obj.Hv(d, Hd);

		alpha = rTr/ddot_(&n, d, &inc, Hd, &inc);
		daxpy_(&n, &alpha, d, &inc, s, &inc);
		if (dnrm2_(&n, s, &inc) > delta)
		{
			info("cg reaches trust region boundary\n");
			alpha = -alpha;
			daxpy_(&n, &alpha, d, &inc, s, &inc);

			double std = ddot_(&n, s, &inc, d, &inc);
			double sts = ddot_(&n, s, &inc, s, &inc);
			double dtd = ddot_(&n, d, &inc, d, &inc);
			double dsq = delta*delta;
			double rad = sqrt(std*std + dtd*(dsq-sts));
			if (std >= 0)
				alpha = (dsq - sts)/(std + rad);
			else
				alpha = (rad - std)/dtd;
			daxpy_(&n, &alpha, d, &inc, s, &inc);
			alpha = -alpha;
			daxpy_(&n, &alpha, Hd, &inc, r, &inc);
			break;
		}
		alpha = -alpha;
		daxpy_(&n, &alpha, Hd, &inc, r, &inc);
		rnewTrnew = ddot_(&n, r, &inc, r, &inc);
		beta = rnewTrnew/rTr;
		dscal_(&n, &beta, d, &inc);
		daxpy_(&n, &one, r, &inc, d, &inc);
		rTr = rnewTrnew;
	}

	return(cg_iter);
}

double TRON::norm_inf(int n, double *x)
{
	double dmax = fabs(x[0]);
	for (int i=1; i<n; i++)
		if (fabs(x[i]) >= dmax)
			dmax = fabs(x[i]);
	return(dmax);
}

void TRON::set_print_string(void (*print_string) (const char *buf))
{
	tron_print_string = print_string;
}
