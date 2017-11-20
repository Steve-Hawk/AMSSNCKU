//$Id: bssn_rhs.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef BSSN_H
#define BSSN_H

#ifdef fortran1
#define f_compute_rhs_bssn compute_rhs_bssn
#define f_compute_rhs_bssn_ss compute_rhs_bssn_ss
#endif
#ifdef fortran2
#define f_compute_rhs_bssn COMPUTE_RHS_BSSN
#define f_compute_rhs_bssn_ss COMPUTE_RHS_BSSN_SS
#endif
#ifdef fortran3
#define f_compute_rhs_bssn compute_rhs_bssn_
#define f_compute_rhs_bssn_ss compute_rhs_bssn_ss_
#endif
extern "C" { int f_compute_rhs_bssn(int *, double &,double *, double *,double *,  //ex,T,X,Y,Z
		                double *, double *,                     //chi, trK
                                double *, double *,double *, double *, double *,double *, //gij
                                double *, double *,double *, double *, double *,double *, //Aij
                                double *, double *,double *,                              //Gam
				double *, double *, double *,double *, double *, double *,double *, //Gauge
		                double *, double *,                     //chi, trK
                                double *, double *,double *, double *, double *,double *, //gij
                                double *, double *,double *, double *, double *,double *, //Aij
                                double *, double *,double *,                              //Gam
				double *, double *, double *,double *, double *, double *,double *, //Gauge
			       	double *, double *,double *, double *, double *, double *,double *, double *, double *,double *, //stress-energy
                                double *, double *,double *, double *, double *,double *, // Christoffel
                                double *, double *,double *, double *, double *,double *, // Christoffel
                                double *, double *,double *, double *, double *,double *, // Christoffel
                                double *, double *,double *, double *, double *,double *, // Ricci
				int &,int &,double &);}

extern "C" { int f_compute_rhs_bssn_ss(int *, double &,double *, double *,double *,  //ex,T,rho,sigma,R
                                double *, double *, double *,  //X,Y,Z
                		double *,double *,double *,    //drhodx,drhody,drhodz
				double *,double *,double *,    //dsigmadx,dsigmady,dsigmadz
				double *,double *,double *,    //dRdx,dRdy,dRdz
		                double *,double *,double *,double *,double *,double *, //drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
                		double *,double *,double *,double *,double *,double *, //dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
		                double *,double *,double *,double *,double *,double *, //dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
		                double *, double *,                     //chi, trK
                                double *, double *,double *, double *, double *,double *, //gij
                                double *, double *,double *, double *, double *,double *, //Aij
                                double *, double *,double *,                              //Gam
				double *, double *, double *,double *, double *, double *,double *, //Gauge
		                double *, double *,                     //chi, trK
                                double *, double *,double *, double *, double *,double *, //gij
                                double *, double *,double *, double *, double *,double *, //Aij
                                double *, double *,double *,                              //Gam
				double *, double *, double *,double *, double *, double *,double *, //Gauge
			       	double *, double *,double *, double *, double *, double *,double *, double *, double *,double *, //stress-energy
                                double *, double *,double *, double *, double *,double *, // Christoffel
                                double *, double *,double *, double *, double *,double *, // Christoffel
                                double *, double *,double *, double *, double *,double *, // Christoffel
                                double *, double *,double *, double *, double *,double *, // Ricci
				int &,int &,double &,int&);}
#endif    /* BSSN_H */

