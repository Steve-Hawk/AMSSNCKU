//$Id: var.h,v 1.1.1.1 2012/02/03 08:46:28 zjcao Exp $
#ifndef VAR_H
#define VAR_H

class var {

public:
   char name[20];
   int sgfn;
   double SoA[3];

public:
   
   var(char *namei, int sgfni,
       const double SYM1, const double SYM2, const double SYM3);

   ~var();
};

#endif   /* VAR_H */
