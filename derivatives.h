// $Id: derivatives.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef DERIVATIVES
#define DERIVATIVES

#ifdef fortran1
#define f_fderivs fderivs
#define f_fdderivs fdderivs
#endif
#ifdef fortran2
#define f_fderivs FDERIVS
#define f_fdderivs FDDERIVS
#endif
#ifdef fortran3
#define f_fderivs fderivs_
#define f_fdderivs fdderivs_
#endif

extern "C" { void f_fderivs(int *, double *, 
		           double *, double *, double *,
			   double *, double *, double *,
                           double &, double &, double &, int &,int &); }

extern "C" { void f_fdderivs(int *, double *, 
		           double *, double *, double *, double *, double *, double *,
			   double *,double *,double *,
                           double &, double &, double &, int &,int &); }

#endif    /* DERIVATIVES */
