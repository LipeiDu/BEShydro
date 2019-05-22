/*
* HydroPlugin.c
*
*  Created on: Oct 23, 2015
*      Author: bazow
*/
#include <stdlib.h>
// for printf
#include <stdio.h>

// for timing
#include <ctime>
#include <iostream>

//for cornelius and writing freezeout file
#include <fstream>
#include "cornelius-c++-1.3/cornelius.cpp"
#include "FreezeOut.cpp"
#include "Memory.cpp"

//for jet-medium interactions
#include "ToyJetClass.cpp"

#include "edu/osu/rhic/harness/hydro/HydroPlugin.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"
#include "edu/osu/rhic/harness/io/FileIO.h"
#include "edu/osu/rhic/trunk/ic/InitialConditions.h"
#include "edu/osu/rhic/trunk/hydro/FullyDiscreteKurganovTadmorScheme.h"
#include "edu/osu/rhic/trunk/hydro/EnergyMomentumTensor.h"
#include "edu/osu/rhic/trunk/eos/EquationOfState.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalSources.h"//lipei

#define FREQ 100 //write output to file every FREQ timesteps
#define FOFREQ 10 //call freezeout surface finder every FOFREQ timesteps
#define FOTEST 0 //if true, freezeout surface file is written with proper times rounded (down) to step size
#define FOFORMAT 0 // 0 : write f.o. surface to ASCII file ;  1 : write to binary file
#define JET 0 // 0 to turn off jet evolution, 1 to turn it on

void outputDynamicalQuantities(double t, const char *outputDir, void * latticeParams)
{
  output(e, t, outputDir, "e", latticeParams);
  output(p, t, outputDir, "p", latticeParams);
  //output(u->ux, t, outputDir, "ux", latticeParams);
  //output(u->uy, t, outputDir, "uy", latticeParams);
  output(u->un, t, outputDir, "un", latticeParams);
  //output(u->ut, t, outputDir, "ut", latticeParams);
  //output(q->ttt, t, outputDir, "ttt", latticeParams);
  //output(q->ttn, t, outputDir, "ttn", latticeParams);
  //output(termX, t, outputDir, "termx", latticeParams);
  //output(termY, t, outputDir, "termy", latticeParams);
  //output(termZ, t, outputDir, "termz", latticeParams);
  #ifdef PIMUNU
  //output(q->pixx, t, outputDir, "pixx", latticeParams);
  //output(q->pixy, t, outputDir, "pixy", latticeParams);
  //output(q->pixn, t, outputDir, "pixn", latticeParams);
  //output(q->piyy, t, outputDir, "piyy", latticeParams);
  //output(q->piyn, t, outputDir, "piyn", latticeParams);
  //output(q->pinn, t, outputDir, "pinn", latticeParams);
  #endif
  #ifdef PI
  //output(q->Pi, t, outputDir, "Pi", latticeParams);
  #endif
  #ifdef NBMU
  output(rhob, t, outputDir, "rhob", latticeParams);
  output(muB, t, outputDir, "muBT", latticeParams);
  //output(T, t, outputDir, "T", latticeParams);
  //output(termX, t, outputDir, "kappaB", latticeParams);
  //output(term2, t, outputDir, "muB", latticeParams);
  //output(q->Nbt, t, outputDir, "Nbt", latticeParams);
  #endif
  #ifdef VMU
  //output(q->nbt, t, outputDir, "nbtau", latticeParams);
  //output(q->nbx, t, outputDir, "nbx", latticeParams);
  //output(q->nby, t, outputDir, "nby", latticeParams);
  output(q->nbn, t, outputDir, "nbn", latticeParams);
  #endif
}

void outputAnalysis(double t, const char *outputDir, void * latticeParams)
{
    FILE *fp;
    char fname[255];
    sprintf(fname, "%s/AnalysisData.dat", outputDir);
    fp=fopen(fname, "a+");
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double x,y,z;
    
    int i,j,k;
    int s;
    
    double v2t,becc,eecc,v2t1,v2t2;
    v2t = 0;
    v2t1 = 0;
    v2t2 = 0;
    becc = 0;
    eecc = 0;
    double bymx = 0;
    double bxy = 0;
    double bypx = 0;
    double eymx = 0;
    double exy = 0;
    double eypx = 0;
    
    //k=(nz+3)/2;
    //j=(ny+3)/2;
    //i=(nx+3)/2;
    for(k = 2; k < nz+2; ++k) {
        z = (k-2 - (nz-1)/2.)*dz;
        for(j = 2; j < ny+2; ++j) {
            y = (j-2 - (ny-1)/2.)*dy;
            for(i = 2; i < nx+2; ++i) {
                x = (i-2 - (nx-1)/2.)*dx;
                s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                //double tt=Ttt(e[s],p[s],u->ut[s],q->pitt[s]);
                //double tx=Ttx(e[s],p[s],u->ut[s],u->ux[s],q->pitx[s]);
                //double ty=Tty(e[s],p[s],u->ut[s],u->uy[s],q->pity[s]);
                //double tn=Ttn(e[s],p[s],u->ut[s],u->un[s],q->pitn[s]);
#ifndef PIMUNU
                double pixx=0;
                double piyy=0;
#else
                double pixx=q->pixx[s];
                double piyy=q->piyy[s];
#endif
                double xx=Txx(e[s],p[s],u->ux[s],pixx);
                //double xy=Txy(e[s],p[s],u->ux[s],u->uy[s],q->pixy[s]);
                //double xn=Txn(e[s],p[s],u->ux[s],u->un[s],q->pixn[s]);
                double yy=Tyy(e[s],p[s],u->uy[s],piyy);
                //double yn=Tyn(e[s],p[s],u->uy[s],u->un[s],q->piyn[s]);
                //double nn=Tnn(e[s],p[s],u->un[s],q->pinn[s],t);
                
                bymx = bymx + (y*y - x*x)*rhob[s];
                bxy = bxy + x*y*rhob[s];
                bypx = bypx + (y*y + x*x)*rhob[s];
                eymx = eymx + (y*y - x*x)*e[s];
                exy = exy + x*y*e[s];
                eypx = eypx + (y*y + x*x)*e[s];
                v2t1 = v2t1 + (xx-yy);
                v2t2 = v2t2 + (xx+yy);
            }
        }
    }
    
    becc = sqrt(bymx*bymx + 4*bxy*bxy)/bypx;
    eecc = sqrt(eymx*eymx + 4*exy*exy)/eypx;
    v2t = v2t1/v2t2;
    
    fprintf(fp, "%.3f\t%.8f\t%.8f\t%.8f\n",t,v2t,becc,eecc);
    
    fclose(fp);
}


void run(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, const char *outputDir)
{
  struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
  struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
  struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

  //************************************************************************************\
  //* System configuration
  //************************************************************************************/
  int nt = lattice->numProperTimePoints;
  int nx = lattice->numLatticePointsX;
  int ny = lattice->numLatticePointsY;
  int nz = lattice->numLatticePointsRapidity;

  int ncx = lattice->numComputationalLatticePointsX;
  int ncy = lattice->numComputationalLatticePointsY;
  int ncz = lattice->numComputationalLatticePointsRapidity;
  int nElements = ncx * ncy * ncz;

  double t0 = hydro->initialProperTimePoint;
  double dt = lattice->latticeSpacingProperTime;
  double dx = lattice->latticeSpacingX;
  double dy = lattice->latticeSpacingY;
  double dz = lattice->latticeSpacingRapidity;
  double e0 = initCond->initialEnergyDensity;
    
  int initialConditionType = initCond->initialConditionType;//Lipei
  int numberOfSourceFiles = initCond->numberOfSourceFiles;//Lipei

  double freezeoutTemperatureGeV = hydro->freezeoutTemperatureGeV;
  const double hbarc = 0.197326938;
  const double freezeoutTemperature = freezeoutTemperatureGeV/hbarc;
  //const double freezeoutEnergyDensity = e0*pow(freezeoutTemperature,4);
  const double freezeoutEnergyDensity = equilibriumEnergyDensity(freezeoutTemperature);
  //const double freezeoutEnergyDensity = 1.8;//Lipei
    
  printf("Grid size = %d x %d x %d\n", nx, ny, nz);
  printf("spatial resolution = (%.3f, %.3f, %.3f)\n", lattice->latticeSpacingX, lattice->latticeSpacingY, lattice->latticeSpacingRapidity);
  printf("freezeout temperature = %.3f [fm^-1] (eF = %.3f [fm^-4])\n", freezeoutTemperature, freezeoutEnergyDensity);

  // allocate memory
  allocateHostMemory(nElements);
  // Read in the table of Equation of State
#ifdef NBMU
  getEquationOfStateTable();//Lipei
  //testEOS();
#endif

  //************************************************************************************\
  //* Jet stuff
  //************************************************************************************/
  jetParton parton;
  if (JET)
  {
    //declare a jet parton instance
    printf("initializing jet parton at center of grid with momentum in y direction \n");
    //initialize four momenta and position to zero
    for (int ip = 0; ip < 4; ip++)
    {
      parton.momentum[ip] = 0.0;
      parton.position[ip] = 0.0;
    }
    // initialize jet at center of coordinate grid with momenta along y direction
    parton.mass = 1;
    parton.position[0] = t0; //same as hydro start time
    parton.momentum[0] = 12.0; //nonzero p^tau
    parton.momentum[2] = 10.0; //nonzero p^y
  }

  //************************************************************************************\
  //* initialize cornelius for freezeout surface finding
  //************************************************************************************/
  //see example_4d() in example_cornelius
  //this works only for full 3+1 d simulation? need to find a way to generalize to n+1 d
  int dim;
  double *lattice_spacing;
  if ((nx > 1) && (ny > 1) && (nz > 1))
  {
    dim = 4;
    lattice_spacing = new double[dim];
    lattice_spacing[0] = dt;
    lattice_spacing[1] = dx;
    lattice_spacing[2] = dy;
    lattice_spacing[3] = dz;
  }
  else if ((nx > 1) && (ny > 1) && (nz < 2))
  {
    dim = 3;
    lattice_spacing = new double[dim];
    lattice_spacing[0] = dt;
    lattice_spacing[1] = dx;
    lattice_spacing[2] = dy;
  }
  else
  {
    printf("simulation is not in 3+1D or 2+1D; freezeout finder will not work!\n");
  }

  Cornelius cor;
  cor.init(dim, freezeoutEnergyDensity, lattice_spacing);

  double ****energy_density_evoution;
  energy_density_evoution = calloc4dArray(energy_density_evoution, FOFREQ+1, nx, ny, nz);

  //make an array to store all the hydrodynamic variables for FOFREQ time steps
  //to be written to file once the freezeout surface is determined by the critical energy density
  int n_hydro_vars = 16; //u0, u1, u2, u3, e, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, Pi, the temperature and pressure are calclated with EoS
  double *****hydrodynamic_evoution;
  hydrodynamic_evoution = calloc5dArray(hydrodynamic_evoution, n_hydro_vars, FOFREQ+1, nx, ny, nz);

  //for 3+1D simulations
  double ****hyperCube4D;
  hyperCube4D = calloc4dArray(hyperCube4D, 2, 2, 2, 2);
  //for 2+1D simulations
  double ***hyperCube3D;
  hyperCube3D = calloc3dArray(hyperCube3D, 2, 2, 2);

  //open the freezeout surface file
  ofstream freezeoutSurfaceFile;
  if (FOFORMAT == 0) freezeoutSurfaceFile.open("output/surface.dat");
  else freezeoutSurfaceFile.open("output/surface.dat", ios::binary);

  //************************************************************************************\
  //* Fluid dynamic initialization
  //************************************************************************************/
  double t = t0;
  // Generate initial conditions
  setInitialConditions(latticeParams, initCondParams, hydroParams, rootDirectory);
  // Calculate conserved quantities
  setConservedVariables(t, latticeParams);
  // Impose boundary conditions with ghost cells
  setGhostCells(q,e,p,u,latticeParams,rhob,muB,T);//rhob by Lipei

  //************************************************************************************\
  //* Evolve the system in time
  //************************************************************************************/
  int ictr = (nx % 2 == 0) ? ncx/2 : (ncx-1)/2;
  int jctr = (ny % 2 == 0) ? ncy/2 : (ncy-1)/2;
  int kctr = (nz % 2 == 0) ? ncz/2 : (ncz-1)/2;
  int sctr = columnMajorLinearIndex(ictr, jctr, kctr, ncx, ncy);

  std::clock_t t1,t2;

  double totalTime = 0;
  int nsteps = 0;

  int accumulator1 = 0;
  int accumulator2 = 0;
    
  // evolve in time
  for (int n = 1; n <= nt+1; ++n)
  {
    // copy variables back to host and write to disk
    if ((n-1) % FREQ == 0)
    {
      printf("n = %d:%d (t = %.3f),\t (e, p) = (%.3f, %.3f) [fm^-4],\t (rhob = %.3f ),\t (T = %.3f [GeV]),\n",
      n - 1, nt, t, e[sctr], p[sctr], rhob[sctr], effectiveTemperature(e[sctr], rhob[sctr]) * hbarc);
      outputDynamicalQuantities(t, outputDir, latticeParams);
      outputAnalysis(t, outputDir, latticeParams);
      // end hydrodynamic simulation if the temperature is below the freezeout temperature
      //if(e[sctr] < freezeoutEnergyDensity) {
      //printf("\nReached freezeout temperature at the center.\n");
      //break;
      //}
    }

    //************************************************************************************\
    // Freeze-out finder (Derek)
    // the freezeout surface file is written in the format which can
    // be read by iS3D : https://github.com/derekeverett/iS3D
    //************************************************************************************/
    //append the energy density and all hydro variables to storage arrays
    int nFO = n % FOFREQ;

    if(nFO == 0) //swap in the old values so that freezeout volume elements have overlap between calls to finder
    {
      swapAndSetHydroVariables(energy_density_evoution, hydrodynamic_evoution, q, e, u, nx, ny, nz, FOFREQ);
    }
    else //update the values of the rest of the array with current time step
    {
      setHydroVariables(energy_density_evoution, hydrodynamic_evoution, q, e, u, nx, ny, nz, FOFREQ, n);
    }

    //the n=1 values are written to the it = 2 index of array, so don't start until here
    int start;
    if (n <= FOFREQ) start = 2;
    else start = 0;
    if (nFO == FOFREQ - 1) //call the freezeout finder should this be put before the values are set?
    {

      //besides writing centroid and normal to file, write all the hydro variables
      int dimZ;
      if (dim == 4) dimZ = nz-1; //enter the loop over iz, and avoid problems at boundary
      else if (dim == 3) dimZ = 1; //we need to enter the 'loop' over iz rather than skipping it
      for (int it = start; it < FOFREQ; it++) //note* avoiding boundary problems (reading outside array)
      {
        for (int ix = 0; ix < nx-1; ix++)
        {
          for (int iy = 0; iy < ny-1; iy++)
          {
            for (int iz = 0; iz < dimZ; iz++)
            {
              //write the values of energy density to all corners of the hyperCube
              if (dim == 4) writeEnergyDensityToHypercube4D(hyperCube4D, energy_density_evoution, it, ix, iy, iz);
              else if (dim == 3) writeEnergyDensityToHypercube3D(hyperCube3D, energy_density_evoution, it, ix, iy);

              //use cornelius to find the centroid and normal vector of each hyperCube
              if (dim == 4) cor.find_surface_4d(hyperCube4D);
              else if (dim == 3) cor.find_surface_3d(hyperCube3D);
              //write centroid and normal of each surface element to file
              for (int i = 0; i < cor.get_Nelements(); i++)
              {
                double temp = 0.0; //temporary variable
                double cell_tau = t0 + ((double)(n - FOFREQ + it)) * dt; //check if this is the correct time!
                double cell_x = (double)ix * dx  - (((double)(nx-1)) / 2.0 * dx);
                double cell_y = (double)iy * dy  - (((double)(ny-1)) / 2.0 * dy);
                double cell_z = (double)iz * dz  - (((double)(nz-1)) / 2.0 * dz);

                double tau_frac = cor.get_centroid_elem(i,0) / lattice_spacing[0];
                double x_frac = cor.get_centroid_elem(i,1) / lattice_spacing[1];
                double y_frac = cor.get_centroid_elem(i,2) / lattice_spacing[2];
                double z_frac;
                if (dim == 4) z_frac = cor.get_centroid_elem(i,3) / lattice_spacing[3];
                else z_frac = 0.0;

                if (FOFORMAT == 0) //write ASCII file
                {
                  //first write the contravariant position vector
                  if (FOTEST) {freezeoutSurfaceFile << cell_tau << " ";}
                  else {freezeoutSurfaceFile << cor.get_centroid_elem(i,0) + cell_tau << " ";}
                  freezeoutSurfaceFile << cor.get_centroid_elem(i,1) + cell_x << " ";
                  freezeoutSurfaceFile << cor.get_centroid_elem(i,2) + cell_y << " ";
                  if (dim == 4) freezeoutSurfaceFile << cor.get_centroid_elem(i,3) + cell_z << " ";
                  else freezeoutSurfaceFile << cell_z << " ";
                  //then the contravariant surface normal element; note jacobian factors of tau for milne coordinates
                  freezeoutSurfaceFile << t * cor.get_normal_elem(i,0) << " ";
                  freezeoutSurfaceFile << t * cor.get_normal_elem(i,1) << " ";
                  freezeoutSurfaceFile << t * cor.get_normal_elem(i,2) << " ";
                  if (dim == 4) freezeoutSurfaceFile << t * cor.get_normal_elem(i,3) << " ";
                  else freezeoutSurfaceFile << 0.0 << " ";
                  //write all the necessary hydro dynamic variables by first performing linear interpolation from values at
                  //corners of hypercube

                  if (dim == 4) // for 3+1D
                  {
                    //first write the contravariant flow velocity
                    for (int ivar = 0; ivar < dim; ivar++)
                    {
                      temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    //write the energy density
                    temp = interpolateVariable4D(hydrodynamic_evoution, 4, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                    freezeoutSurfaceFile << temp << " "; //note : iSpectra reads in file in fm^x units e.g. energy density should be written in fm^-4
                    //the temperature !this needs to be checked
                    freezeoutSurfaceFile << effectiveTemperature(temp) << " ";
                    //the thermal pressure
                    freezeoutSurfaceFile << equilibriumPressure(temp) << " ";
                    //write ten components of pi_(mu,nu) shear viscous tensor
                    for (int ivar = 5; ivar < 15; ivar++)
                    {
                      temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    //write the bulk pressure Pi, and start a new line
                    temp = interpolateVariable4D(hydrodynamic_evoution, 15, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                    freezeoutSurfaceFile << temp << endl;
                  }

                  else //for 2+1D
                  {
                    //first write the contravariant flow velocity
                    for (int ivar = 0; ivar < 4; ivar++)
                    {
                      temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    //write the energy density
                    temp = interpolateVariable3D(hydrodynamic_evoution, 4, it, ix, iy, tau_frac, x_frac, y_frac);
                    freezeoutSurfaceFile << temp << " "; //note units of fm^-4 appropriate for iSpectra reading
                    //the temperature !this needs to be checked
                    freezeoutSurfaceFile << effectiveTemperature(temp) << " ";
                    //the thermal pressure
                    freezeoutSurfaceFile << equilibriumPressure(temp) << " ";
                    //write ten components of pi_(mu,nu) shear viscous tensor
                    for (int ivar = 5; ivar < 15; ivar++)
                    {
                      temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    //write the bulk pressure Pi, and start a new line
                    temp = interpolateVariable3D(hydrodynamic_evoution, 15, it, ix, iy, tau_frac, x_frac, y_frac);
                    freezeoutSurfaceFile << temp << endl;
                  }
                }

                /*
                else //write in binary
                {
                  if (FOTEST) {freezeoutSurfaceFile.write(cell_tau, sizeof(double));}
                  else {freezeoutSurfaceFile.write(cor.get_centroid_elem(i,0) + cell_tau, sizeof(double));}
                  freezeoutSurfaceFile.write(cor.get_centroid_elem(i,1) + cell_x,sizeof(double));
                  freezeoutSurfaceFile.write(cor.get_centroid_elem(i,2) + cell_y,sizeof(double));
                  if (dim == 4) freezeoutSurfaceFile.write(cor.get_centroid_elem(i,3) + cell_z,sizeof(double));
                  else freezeoutSurfaceFile.write(cell_z,sizeof(double));
                  //then the surface normal element; note jacobian factors of +/- tau for milne coordinates
                  freezeoutSurfaceFile.write(t * cor.get_normal_elem(i,0),sizeof(double));
                  freezeoutSurfaceFile.write((-1.0) * t * cor.get_normal_elem(i,1),sizeof(double));
                  freezeoutSurfaceFile.write((-1.0) * t * cor.get_normal_elem(i,2),sizeof(double));
                  if (dim == 4) freezeoutSurfaceFile.write((-1.0) * t * cor.get_normal_elem(i,3),sizeof(double));
                  else freezeoutSurfaceFile.write(0.0,sizeof(double));
                  //write all the necessary hydro dynamic variables by first performing linear interpolation from values at
                  //corners of hypercube

                  if (dim == 4) // for 3+1D
                  {
                    //first write the flow velocity
                    for (int ivar = 0; ivar < dim; ivar++)
                    {
                      temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                      freezeoutSurfaceFile.write(temp, sizeof(double));
                    }
                    //write the energy density
                    temp = interpolateVariable4D(hydrodynamic_evoution, 4, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                    freezeoutSurfaceFile.write(temp, sizeof(double));
                    //the temperature !this needs to be checked
                    freezeoutSurfaceFile.write(effectiveTemperature(temp), sizeof(double));
                    //the baryon chemical potential, writing zero for now
                    freezeoutSurfaceFile.write(0.0, sizeof(double));
                    //  (e + P) / T , the entropy density for zero chem. potentials !check this, note we could be a divide by zero problem if T=0!
                    double e_plus_P_over_T = (temp + equilibriumPressure(temp)) / effectiveTemperature(temp);
                    freezeoutSurfaceFile.write(e_plus_P_over_T, sizeof(double));
                    //write ten components of pi_(mu,nu) shear viscous tensor
                    for (int ivar = 5; ivar < 15; ivar++)
                    {
                      temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                      freezeoutSurfaceFile.write(temp, sizeof(double));
                    }
                    //write the bulk pressure Pi, and start a new line
                    temp = interpolateVariable4D(hydrodynamic_evoution, 15, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                    freezeoutSurfaceFile.write(temp, sizeof(double));
                  }

                  else //for 2+1D
                  {
                    //first write the flow velocity
                    for (int ivar = 0; ivar < 4; ivar++)
                    {
                      temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                      freezeoutSurfaceFile.write(temp, sizeof(double));
                    }
                    //write the energy density
                    temp = interpolateVariable3D(hydrodynamic_evoution, 4, it, ix, iy, tau_frac, x_frac, y_frac);
                    freezeoutSurfaceFile.write(temp, sizeof(double)); //note factors of hbarc to give units (GeV/fm^3)
                    //the temperature !this needs to be checked
                    freezeoutSurfaceFile.write(effectiveTemperature(temp), sizeof(double));
                    //the baryon chemical potential, writing zero for now
                    freezeoutSurfaceFile.write(0.0, sizeof(double));
                    //  (e + P) / T , the entropy density for zero chem. potentials !check this, note we could be a divide by zero problem if T=0!
                    double e_plus_P_over_T = (temp + equilibriumPressure(temp)) / effectiveTemperature(temp);
                    freezeoutSurfaceFile.write(e_plus_P_over_T, sizeof(double));
                    //write ten components of pi_(mu,nu) shear viscous tensor
                    for (int ivar = 5; ivar < 15; ivar++)
                    {
                      temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                      freezeoutSurfaceFile.write(temp, sizeof(double));
                    }
                    //write the bulk pressure Pi, and start a new line
                    temp = interpolateVariable3D(hydrodynamic_evoution, 15, it, ix, iy, tau_frac, x_frac, y_frac);
                    freezeoutSurfaceFile.write(temp, sizeof(double));
                  }
                }
                */
              }
            }
          }
        }
      }
    }
    
    //if all cells are below freezeout temperature end hydro
    accumulator1 = 0;
    for (int ix = 2; ix < nx+2; ix++)
    {
      for (int iy = 2; iy < ny+2; iy++)
      {
        for (int iz = 2; iz < nz+2; iz++)
        {
          int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4);
          if (e[s] > freezeoutEnergyDensity) accumulator1 = accumulator1 + 1;
        }
      }
    }
    if (accumulator1 == 0) accumulator2 += 1;
    if (accumulator2 >= FOFREQ+1 && (n > numberOfSourceFiles)) //only break once freezeout finder has had a chance to search/write to file
    {
      printf("\nAll cells have dropped below freezeout energy density\n");
      break;
    }
    
    //************************************************************************************\
    //* Evolution by RungeKutta
    //************************************************************************************/

    t1 = std::clock();


    //****************JET STUFF************/
    if (JET)
    {
      //get the local fluid velocity and energy density/temperature and evolve jet momentum
      parton.energyLoss(nx, ny, nz, dt, dx, dy, dz, u->ut, u->ux, u->uy, u->un, e, rhob);
      //evolve the jet parton position
      parton.updatePosition(dt);

      //set hydro source terms
      setDynamicalSources(latticeParams, initCondParams, parton.dp_dtau, parton.position);
    }
    //****************JET STUFF************/

    // Read in source terms from particles
    if(initialConditionType==13){
          if(n<=numberOfSourceFiles) readInSource(n, latticeParams, initCondParams, hydroParams, rootDirectory);
          else noSource(latticeParams, initCondParams);
    }

    rungeKutta2(t, dt, q, Q, latticeParams, hydroParams);
    t2 = std::clock();
    double delta_time = (t2 - t1) / (double)(CLOCKS_PER_SEC / 1000);
    if ((n-1) % FREQ == 0) printf("(Elapsed time: %.3f ms)\n",delta_time);
    totalTime+=delta_time;
    ++nsteps;

    setCurrentConservedVariables();

    t = t0 + n * dt;
  }
  printf("Average time/step: %.3f ms\n",totalTime/((double)nsteps));

  freezeoutSurfaceFile.close();
  //************************************************************************************\
  //* Deallocate host memory
  //************************************************************************************/
  freeHostMemory();

  //Deallocate memory used for freezeout finding
  free4dArray(energy_density_evoution, FOFREQ+1, nx, ny);
  free5dArray(hydrodynamic_evoution, n_hydro_vars, FOFREQ+1, nx, ny);
  delete [] lattice_spacing;

  free4dArray(hyperCube4D, 2, 2, 2);
}
