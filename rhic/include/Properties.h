//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#include <libconfig.h>

void getIntegerProperty(config_t *cfg, const char* propName, int *propValue, int defaultValue);
void getDoubleProperty(config_t *cfg, const char* propName, double *propValue, double defaultValue);

#endif /* PROPERTIES_H_ */
