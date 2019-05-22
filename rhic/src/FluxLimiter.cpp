/*
 * FluxLimiter.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

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
    /*PRECISION l = THETA * (y - x);
    PRECISION c = (z - x) / 2;
    PRECISION r = THETA * (z - y);
    return minmod3(l, c, r);*/
    
    //Lipei
    //MUSIC
    PRECISION l = THETA * (y - x);
    PRECISION c = (z - x) / 2;
    PRECISION r = THETA * (z - y);
    double tem;
    if( (l > 0.0) && (c > 0.0) && (r > 0.0) )
    {
        tem = mini(l, c);
        return mini(r, tem);
    }
    else if( (l < 0.0) && (c < 0.0) && (r < 0.0) )
    {
        tem = maxi(l, c);
        return maxi(r, tem);
    }
    else
        return 0.0;
    
    //UMIST
    /*PRECISION min1 = mini(2*(y-x), (0.25*(z-y)+0.75*(y-x)));
    PRECISION min2 = mini((0.75*(z-y)+0.25*(y-x)),2*(z-y));
    PRECISION min3 = mini(min1,min2);
    return maxi(0, min3);*/
    
    //van Albada2
    //PRECISION pr=(y-x)/(z-y);
    //return 2*(y-x)*(z-y)/((y-x)*(y-x)+(z-y)*(z-y));
    
    //ospre
    /*PRECISION pr=(y-x)/(z-y);
    return 1.5*(pr*pr+pr)/(pr*pr+pr+1);*/
    
    //superbee
    /*PRECISION min1 = mini(2*(y-x),(z-y));
    PRECISION min2 = mini(y-x, 2*(z-y));
    PRECISION max1 = maxi(0, min1);
    return maxi(max1, min2);*/
    
    //smart
    /*PRECISION l = 2*(y-x);
    PRECISION c = 0.25*(z-y)+0.75*(y-x);
    PRECISION r = 4*(z-y);
    double tem;
    if( (l > 0.0) && (c > 0.0) && (r > 0.0) )
    {
        tem = mini(l, c);
        return mini(r, tem);
    }
    else if( (l < 0.0) && (c < 0.0) && (r < 0.0) )
    {
        tem = maxi(l, c);
        return maxi(r, tem);
    }
    else
        return 0.0;*/
    
    //Koren
    /*PRECISION l = 2*(y-x);
    PRECISION c = 0.333333*((z-y)+2*(y-x));
    PRECISION r = 2*(z-y);
    double tem;
    if( (l > 0.0) && (c > 0.0) && (r > 0.0) )
    {
        tem = mini(l, c);
        return mini(r, tem);
    }
    else if( (l < 0.0) && (c < 0.0) && (r < 0.0) )
    {
        tem = maxi(l, c);
        return maxi(r, tem);
    }
    else
        return 0.0;*/
    
    //monotonized central
    /*PRECISION l = 2*(y-x);
    PRECISION c = 0.5*((z-y)+(y-x));
    PRECISION r = 2*(z-y);
    double tem;
    if( (l > 0.0) && (c > 0.0) && (r > 0.0) )
    {
        tem = mini(l, c);
        return mini(r, tem);
    }
    else if( (l < 0.0) && (c < 0.0) && (r < 0.0) )
    {
        tem = maxi(l, c);
        return maxi(r, tem);
    }
    else
        return 0.0;*/
    
}
