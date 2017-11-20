//$Id: fmisc.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef FMISC_H
#define FMISC_H

#ifdef fortran1
#define f_pointcopy pointcopy
#define f_copy copy
#define f_global_interp global_interp
#define f_global_interp_ss global_interp_ss
#define f_global_interpind global_interpind
#define f_global_interpind2d global_interpind2d
#define f_global_interpind1d global_interpind1d
#define f_l2normhelper l2normhelper
#define f_average average
#define f_average3 average3
#define f_average2 average2
#define f_average2p average2p
#define f_average2m average2m
#define f_lowerboundset lowerboundset
#endif
#ifdef fortran2
#define f_pointcopy POINTCOPY
#define f_copy COPY
#define f_global_interp GLOBAL_INTERP
#define f_global_interp_ss GLOBAL_INTERP_SS
#define f_global_interpind GLOBAL_INTERPIND
#define f_global_interpind2d GLOBAL_INTERPIND2D
#define f_global_interpind1d GLOBAL_INTERPIND1D
#define f_l2normhelper L2NORMHELPER
#define f_average AVERAGE
#define f_average3 AVERAGE3
#define f_average2 AVERAGE2
#define f_average2p AVERAGE2P
#define f_average2m AVERAGE2M
#define f_lowerboundset LOWERBOUNDSET
#endif
#ifdef fortran3
#define f_pointcopy pointcopy_
#define f_copy copy_
#define f_global_interp global_interp_
#define f_global_interp_ss global_interp_ss_
#define f_global_interpind global_interpind_
#define f_global_interpind2d global_interpind2d_
#define f_global_interpind1d global_interpind1d_
#define f_l2normhelper l2normhelper_
#define f_average average_
#define f_average3 average3_
#define f_average2 average2_
#define f_average2p average2p_
#define f_average2m average2m_
#define f_lowerboundset lowerboundset_
#endif

extern "C" { void f_pointcopy( int &,
		          double*, double *,int*, double*, 
 		          double&,double&,double&,double&); }

extern "C" { void f_copy( int &,
		          double*, double *,int*, double*, 
 		          double*, double *,int*, double*,
			  double*, double*); }

extern "C" { void f_global_interp( int *, double *, double*,double*,
	                           double *, double &,
	                           double &,double &,double &,
			           int &, double *,int &);}

extern "C" { void f_global_interp_ss( int *, double *, double*,double*,
	                           double *, double &,
	                           double &,double &,double &,
			           int &, double *,int &,int &);}

extern "C" { void f_global_interpind( int *, double *, double*,double*,
	                           double *, double &,
	                           double &,double &,double &,
			           int &, double *,int &,
				   int *,double *,int &);}

extern "C" { void f_global_interpind2d( int *, double *, double*,double*,
	                           double *, double &,
	                           double &,double &,double &,
			           int &, double *,int &,
				   int *,double *,int &);}

extern "C" { void f_global_interpind1d( int *, double *, double*,double*,
	                           double *, double &,
	                           double &,double &,double &,
			           int &, double *,int &,
				   int *,double *,int &,int &);}

extern "C" { void f_l2normhelper( int *, double *, double*,double*,
	                           double &, double &,double &,
	                           double &,double &,double &,
			           double *,double &,int &);}

extern "C" { void f_average( int*, double*, double*,double*); }

extern "C" { void f_average3( int*, double*, double*,double*); }

extern "C" { void f_average2( int*, double*, double*, double*,double*); }

extern "C" { void f_average2p( int*, double*, double*, double*,double*); }

extern "C" { void f_average2m( int*, double*, double*, double*,double*); }

extern "C" { void f_lowerboundset(int *, double *,double &);}

#endif    /* FMISC_H */
