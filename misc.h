//$Id: misc.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef MISC_H
#define MISC_H

#ifdef newc
#include <algorithm>   
#include <functional> 
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif

namespace misc
{
inline string&  lTrim(string   &ss)  {   
    string::iterator  p=find_if(ss.begin(),ss.end(),not1(ptr_fun<int,int>(isspace)));   
    ss.erase(ss.begin(),p);   
    return  ss;   
}   
inline string&  rTrim(string   &ss)  {    
    string::reverse_iterator  p=find_if(ss.rbegin(),ss.rend(),not1(ptr_fun<int,int>(isspace)));   
    ss.erase(p.base(),ss.end());   
    return   ss;   
}   
inline string& Trim(string   &st)  {   
    lTrim(rTrim(st));   
    return   st;   
} 

template <typename T>
void swap(T &a,T &b)
{
  T c=a;
  a=b;
  b=c;
}
void tillherecheck(int myrank);
void tillherecheck(const char str[]);
int parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind);
int parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind1, int& ind2);
int parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind1, int& ind2, int& ind3);
void gaulegf(double x1, double x2, double *x, double *w, int n);
void inversearray(double *aa,int NN);
double fact(int N);
double Wigner_d_function(int l,int m, int s, double costheta);
}
#endif   /* MISC_H */
