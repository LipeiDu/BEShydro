//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//


#ifndef FREEZEOUT_H
#define FREEZEOUT_H

#include <stdlib.h>
#include <stdio.h> // for printf

#include "DynamicalVariables.h"

using namespace std;

class FO_Element {
private:
    
public:
    //covariant freezeout cell coordinates
    PRECISION tau;
    PRECISION x;
    PRECISION y;
    PRECISION eta;
    
    //covariant freezeout area normal vector
    PRECISION dat;
    PRECISION dax;
    PRECISION day;
    PRECISION dan;
    
    //contravariant flow velocity
    PRECISION ux;
    PRECISION uy;
    PRECISION un;
    
    //energy density in fm^-4
    PRECISION E;
    //temperature in fm^-1
    PRECISION T;
    //thermal pressure in fm^-4
    PRECISION P;
    
    //five (indep) contravariant components of shear stress
    PRECISION pixx;
    PRECISION pixy;
    PRECISION pixn;
    PRECISION piyy;
    PRECISION piyn;
    
    //the bulk viscous pressure
    PRECISION Pi;
    
#ifdef THERMAL_VORTICITY
    //six components of thermal vorticity tensor
    PRECISION wtx;
    PRECISION wty;
    PRECISION wtn;
    PRECISION wxy;
    PRECISION wxn;
    PRECISION wyn;
#endif
};

#endif
