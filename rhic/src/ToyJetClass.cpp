//#include "ToyJetClass.h"
#include <stdio.h>
#include <math.h>
#include "../include/EquationOfState.h"
#include "../include/DynamicalVariables.h"

class jetParton{
public:
  //contravariant momentum four vector of parton p^{\mu} in milne coords
  double momentum[4];
  //contravariant position four vector of parton x^{\mu} in milne coords
  double position[4];
  //proper time derivative of four momentum
  double dp_dtau[4];
  //mass of parton
  double mass;

  //update the four momentum of parton; just an euler step
  void updateMomentum(double dtau);
  //update the four position of the parton - NOT CORRECT FIX IT
  void updatePosition(double dtau);
  //get parton energy loss based on fluid variables
  void energyLoss(int nx, int ny, int nz, double dt, double dx, double dy, double dz, double *ut, double* ux, double *uy, double* un, double *e, double *rhob);
};

//update the four momentum of parton; just an euler step
void jetParton::updateMomentum(double dtau)
{
  for (int i = 0; i < 4; i++)
  {
    momentum[i] = momentum[i] + (dtau * dp_dtau[i]);
  }
}
//update the four position of the parton - NOT CORRECT FIX IT
void jetParton::updatePosition(double dtau)
{
  position[0] = position[0] + dtau;
  for (int i = 1; i < 4; i++)
  {
    position[i] = position[i] + (dtau * momentum[i] / mass);
  }
}


void jetParton::energyLoss(int nx, int ny, int nz, double dt, double dx, double dy, double dz, double *ut, double* ux, double *uy, double* un, double *e, double *rhob)
{

  int ncx = nx + 4;
  int ncy = ny + 4;

  //get the partons coordinate index
  double xmin = (-1.0) * ((double)(nx - 1) / 2.0) * dx;
  int ix = (int)round((position[1] - xmin) / dx);
  //printf("ix = %d \n", ix);
  double ymin = (-1.0) * ((double)(ny - 1) / 2.0) * dy;
  int iy = (int)round((position[2] - ymin) / dy);
  //printf("iy = %d \n", iy);
  double zmin = (-1.0) * ((double)(nz - 1) / 2.0) * dz;
  int iz = (int)round((position[3] - zmin) / dz);
  //printf("iz = %d \n", iz);
  int s = columnMajorLinearIndex(ix, iy, iz, ncx, ncy);

  //safegaurd in case parton leaves grid!
  if ((ix >= 0) && (ix < nx) && (iy >= 0) && (iy < ny) && (iz >= 0) && (iz < nz))
  {
    //this prescription for drag/energy loss should work for the bjorken flow, may not make sense in general!
    //define a temperature dependent jet-medium 'relaxation time'
    double t_R = 3.0 / (effectiveTemperature(e[s],rhob[s]) * effectiveTemperature(e[s],rhob[s]) * effectiveTemperature(e[s],rhob[s]));

    dp_dtau[0] = (-1.0) * ((momentum[0] / mass) - ut[s]) / t_R;
    dp_dtau[1] = (-1.0) * ((momentum[1] / mass) - ux[s]) / t_R;
    dp_dtau[2] = (-1.0) * ((momentum[2] / mass) - uy[s]) / t_R;
    dp_dtau[3] = (-1.0) * ((momentum[3] / mass) - un[s]) / t_R;

    updateMomentum(dt);
  }

  else
  {
    printf("parton has escaped medium with final position {%f, %f, %f, %f} and momentum {%f, %f, %f, %f} \n", position[0],position[1],position[2],position[3],momentum[0],momentum[1],momentum[2],momentum[3]);
  }

}
