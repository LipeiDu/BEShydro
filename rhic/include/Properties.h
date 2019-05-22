/*
 * Properties.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#include <libconfig.h>

void getIntegerProperty(config_t *cfg, const char* propName, int *propValue, int defaultValue);
void getDoubleProperty(config_t *cfg, const char* propName, double *propValue, double defaultValue);

#endif /* PROPERTIES_H_ */
