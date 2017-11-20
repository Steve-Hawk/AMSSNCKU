//$Id: rungekutta4_rout.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H

#ifdef fortran1
#define f_rungekutta4_rout rungekutta4_rout
#define f_rungekutta4_scalar rungekutta4_scalar
#endif
#ifdef fortran2
#define f_rungekutta4_rout RUNGEKUTTA4_ROUT
#define f_rungekutta4_scalar RUNGEKUTTA4_SCALAR
#endif
#ifdef fortran3
#define f_rungekutta4_rout rungekutta4_rout_
#define f_rungekutta4_scalar rungekutta4_scalar_
#endif

extern "C" { void f_rungekutta4_scalar( double &, double &, double &, double &, int &);}

extern "C" { int  f_rungekutta4_rout( int *, double &, 
			      double*, double *,double *,
			      int &); }

#endif    /* RUNGEKUTTA4_H */
