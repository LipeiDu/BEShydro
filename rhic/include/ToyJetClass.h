//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

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
  void energyLoss(int nx, int ny, int nz, double t, double dt, double dx, double dy, double dz, double *ut, double* ux, double *uy, double* un, double *e, double *rhob, int energyLossType, double initialPartionPositionX);
};
