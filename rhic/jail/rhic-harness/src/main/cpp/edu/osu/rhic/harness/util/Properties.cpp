/*
 * Properties.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/harness/util/Properties.h"

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
