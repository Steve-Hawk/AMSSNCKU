//$Id: getnp4.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef GETNP4_H
#define GETNP4_H

#ifdef fortran1
#define f_getnp4 getnp4
#define f_getnp4_ss getnp4_ss
#endif
#ifdef fortran2
#define f_getnp4 GETNP4
#define f_getnp4_ss GETNP4_SS
#endif
#ifdef fortran3
#define f_getnp4 getnp4_
#define f_getnp4_ss getnp4_ss_
#endif

extern "C" { void f_getnp4(int *, double *, double *,double *,
		                double *, double *,
                                double *, double *,double *, double *, double *,double *,
                                double *, double *,double *, double *, double *,double *,
				double *, double *, double *,double *, double *,double *,
                                double *, double *,double *, double *, double *,double *,
                                double *, double *,double *, double *, double *,double *,
                                double *, double *,double *, double *, double *,double *,
                                double *, double *,int &);}

extern "C" { void f_getnp4_ss(int *, double *, double *,double *,double *, double *,double *,
                                double *, double *,double *,
                                double *, double *,double *,
                                double *, double *,double *,
                                double *, double *,double *,double *, double *,double *,
                                double *, double *,double *,double *, double *,double *,
                                double *, double *,double *,double *, double *,double *,
		                double *, double *,
                                double *, double *,double *, double *, double *,double *,
                                double *, double *,double *, double *, double *,double *,
				double *, double *, double *,double *, double *,double *,
                                double *, double *,double *, double *, double *,double *,
                                double *, double *,double *, double *, double *,double *,
                                double *, double *,double *, double *, double *,double *,
                                double *, double *,int &,int &);}

#endif    /* GETNP4_H */
