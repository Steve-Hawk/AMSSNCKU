//$Id: sommerfeld_rout.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef SOMMERFELD_ROUT_H
#define SOMMERFELD_ROUT_H

#ifdef fortran1
#define f_sommerfeld_rout sommerfeld_rout
#define f_sommerfeld_routbam sommerfeld_routbam
#define f_sommerfeld_routbam_ss sommerfeld_routbam_ss
#endif
#ifdef fortran2
#define f_sommerfeld_rout SOMMERFELD_ROUT
#define f_sommerfeld_rout SOMMERFELD_ROUTBAM
#define f_sommerfeld_rout_ss SOMMERFELD_ROUTBAM_SS
#endif
#ifdef fortran3
#define f_sommerfeld_rout sommerfeld_rout_
#define f_sommerfeld_routbam sommerfeld_routbam_
#define f_sommerfeld_routbam_ss sommerfeld_routbam_ss_
#endif

extern "C" { void f_sommerfeld_rout( int *, double *,double *,double *,
		            double &,double &,double &,double &,double &,double &,double &, double *,
	                    double *,double *,double *,double *,
	                    int &, int &
			    );}

extern "C" { void f_sommerfeld_routbam( int *, double *,double *,double *,
		            double &,double &,double &,double &,double &,double &, double *,
	                    double *,double &,double *, int &
			    );}

extern "C" { void f_sommerfeld_routbam_ss( int *, double *,double *,double *,
		            double &,double &,double &,double &,double &,double &, double *,
	                    double *,double &,double *, int &
			    );}

#endif    /* SOMMERFELD_ROUT_H */
