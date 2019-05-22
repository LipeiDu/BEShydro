/*
 * SourcePart.h
 *
 *  Created on: Nov 25, 2017
 *      Author: Lipei
 */

#ifndef SOURCE_H_
#define SOURCE_H_

void noSource(void * latticeParams, void * initCondParams);
void readInSource(int n, void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory);
void setSource(int n, void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory);
void setDynamicalSources(void * latticeParams, void * initCondParams, double *dp_dtau, double *pos); //function to set dynamical source terms for stress tensor and baryon current
//PRECISION updateDynamicalSources(); //this function reads in source terms from external code
#endif /* SOURCEPART_H_ */
