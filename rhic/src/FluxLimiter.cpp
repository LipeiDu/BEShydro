//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <math.h>

#include "../include/FluxLimiter.h"
#include "../include/DynamicalVariables.h"

#define maxi(a, b) ((a) > (b) ? (a) : (b))
#define mini(a, b) ((a) < (b) ? (a) : (b))

inline int sign(PRECISION x) {
	if (x<0) return -1;
	else return 1;
}
 
inline PRECISION minmod(PRECISION x, PRECISION y) {
	return (sign(x)+sign(y))*fmin(fabs(x),fabs(y))/2;
}

PRECISION minmod3(PRECISION x, PRECISION y, PRECISION z) {
   return minmod(x,minmod(y,z));
}
 
PRECISION approximateDerivative(PRECISION x, PRECISION y, PRECISION z) {
    
    //Generalized minmod limiter
    PRECISION l = THETA * (y - x);
    PRECISION c = (z - x) / 2;
    PRECISION r = THETA * (z - y);
    return minmod3(l, c, r);
}
