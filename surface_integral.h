//$Id: surface_integral.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef SURFACE_INTEGRAL_H
#define SURFACE_INTEGRAL_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <strstream>
#include <math.h>
#endif

#include "cgh.h"
#include "ShellPatch.h"
#include "var.h"
#include "monitor.h"

class surface_integral {

 private: 
   int Symmetry;
  int N_theta, N_phi;           // Number of points in Theta & Phi directions
  double dphi,dcostheta;  
  double *arcostheta, *wtcostheta; 
  int n_tot;                    // size of arrays

  double *nx_g, *ny_g, *nz_g;   // global list of unit normals
  int myrank,cpusize;

 public:
  surface_integral(int iSymmetry);
  ~surface_integral();

  void surf_Wave(double rex,int lev,cgh *GH, var *Rpsi4, var *Ipsi4,
		 int spinw,int maxl, int NN, double *RP, double *IP,
		 monitor *Monitor);  //NN is the length of RP and IP
  void surf_Wave(double rex,int lev,ShellPatch *GH, var *Rpsi4, var *Ipsi4,
		                   int spinw,int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
};
#endif  /* SURFACE_INTEGRAL_H */
