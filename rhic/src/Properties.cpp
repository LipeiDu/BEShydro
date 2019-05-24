//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include "../include/Properties.h"

void getIntegerProperty(config_t *cfg, const char* propName, int *propValue, int defaultValue) {
	  if(config_lookup_int(cfg, propName, propValue))
	    return;
	  else
	    *propValue = defaultValue;
}

void getDoubleProperty(config_t *cfg, const char* propName, double *propValue, double defaultValue) {
	  if(config_lookup_float(cfg, propName, propValue))
	    return;
	  else
	    *propValue = defaultValue;
}
