//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <iostream>
#include <vector>

//for cornelius and writing freezeout file
#include <fstream>
#include "FreezeOut.cpp"
#include "Memory.cpp"

#include "../include/HydroPlugin.h"
#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"
#include "../include/FileIO.h"
#include "../include/InitialConditions.h"
#include "../include/FullyDiscreteKurganovTadmorScheme.h"
#include "../include/PrimaryVariables.h"
#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"

#include "../include/ToyJetClass.h"
#include "../include/HydroPlus.h"
#include "../include/DynamicalSources.h"
#include "../include/HydroAnalysis.h"

#define FREQ 200 //write output to file every FREQ timesteps
#define FOFREQ 10 //call freezeout surface finder every FOFREQ timesteps
#define FOTEST 0 //if true, freezeout surface file is written with proper times rounded (down) to step size
#define JET 0 // 0 to turn off jet evolution, 1 to turn it on

/**************************************************************************************************************************************************/
/* choose dynamical quantities to output
/**************************************************************************************************************************************************/
void outputDynamicalQuantities(double t, const char *outputDir, void * latticeParams)
{
  output(e, t, outputDir, "e", latticeParams);
  output(p, t, outputDir, "p", latticeParams);
  //output(seq, t, outputDir, "seq", latticeParams);
  //output(u->ux, t, outputDir, "ux", latticeParams);
  //output(u->uy, t, outputDir, "uy", latticeParams);
  output(u->un, t, outputDir, "un", latticeParams);
  output(u->ut, t, outputDir, "ut", latticeParams);
  //output(q->ttt, t, outputDir, "ttt", latticeParams);
  //output(q->ttx, t, outputDir, "ttx", latticeParams);
  //output(q->tty, t, outputDir, "tty", latticeParams);
  //output(q->ttn, t, outputDir, "ttn", latticeParams);
  #ifdef PIMUNU
  //output(q->pitx, t, outputDir, "pitx", latticeParams);
  //output(q->pixx, t, outputDir, "pixx", latticeParams);
  //output(q->pixy, t, outputDir, "pixy", latticeParams);
  //output(q->pixn, t, outputDir, "pixn", latticeParams);
  //output(q->piyy, t, outputDir, "piyy", latticeParams);
  //output(q->piyn, t, outputDir, "piyn", latticeParams);
  //output(q->pinn, t, outputDir, "pinn", latticeParams);
  #endif
  #ifdef PI
  output(q->Pi, t, outputDir, "Pi", latticeParams);
  #endif
  output(T, t, outputDir, "T", latticeParams);
  #ifdef NBMU
  output(rhob, t, outputDir, "rhob", latticeParams);
  output(alphaB, t, outputDir, "alphaB", latticeParams);
  //output(q->Nbt, t, outputDir, "Nbt", latticeParams);
  #endif
  #ifdef VMU
  output(q->nbt, t, outputDir, "nbtau", latticeParams);
  //output(q->nbx, t, outputDir, "nbx", latticeParams);
  //output(q->nby, t, outputDir, "nby", latticeParams);
  output(q->nbn, t, outputDir, "nbn", latticeParams);
  //outputPhaseDiagram(alphaB, T, t, outputDir, "muBT", latticeParams);
  #endif
  #ifdef HydroPlus
  //output(q->phiQ[0], t, outputDir, "phiQ0", latticeParams);
  //output(q->phiQ[1], t, outputDir, "phiQ1", latticeParams);
  //output(q->phiQ[2], t, outputDir, "phiQ2", latticeParams);
  //output(eqPhiQ->phiQ[0], t, outputDir, "eqPhiQ0", latticeParams);
  //output(eqPhiQ->phiQ[1], t, outputDir, "eqPhiQ1", latticeParams);
  //output(eqPhiQ->phiQ[2], t, outputDir, "eqPhiQ2", latticeParams);
  #endif
}

/**************************************************************************************************************************************************/
/* the main structure of the hydrodynamic code
/**************************************************************************************************************************************************/
void run(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, const char *outputDir)
{
    
  //************************************************************************************\
  //* System configuration
  //************************************************************************************/
    
  struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
  struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
  struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

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
    
  int initialConditionType = initCond->initialConditionType;
  int numberOfSourceFiles = initCond->numberOfSourceFiles;

  double freezeoutEnergyDensityGeV = hydro->freezeoutEnergyDensityGeV;
  const double hbarc = 0.197326938;
  const double freezeoutEnergyDensity = freezeoutEnergyDensityGeV/hbarc;
    
  printf("Grid size = %d x %d x %d\n", nx, ny, nz);
  printf("spatial resolution = (%.3f, %.3f, %.3f)\n", lattice->latticeSpacingX, lattice->latticeSpacingY, lattice->latticeSpacingRapidity);
  printf("freezeout energy density eF = %.3f [fm^-4]\n", freezeoutEnergyDensity);

  //************************************************************************************\
  //* Allocation and reading in tables
  //************************************************************************************/
    
  // allocate memory
  allocateHostMemory(nElements);
  // Read in the table of Equation of State with baryon
  //getEquationOfStateTable();
  getEquationOfStateTableNEOS();
  //testEOS();
  // baryon diffusion coefficients table
  getBaryonDiffusionCoefficientTable();
  // read in the parameterized correlation length xi(T, muB)
  getCorrelationLengthTable();
  //testCorreLength();

  //************************************************************************************\
  //* Jet initialization
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
    
  std::vector<FO_Element> fo_surf;

  double ****energy_density_evoution;
  energy_density_evoution = calloc4dArray(energy_density_evoution, FOFREQ+1, nx, ny, nz);

  //make an array to store all the hydrodynamic variables for FOFREQ time steps
  //to be written to file once the freezeout surface is determined by the critical energy density
  int n_hydro_vars = 21; //u0, u1, u2, u3, e, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, Pi, rhob, nb0, nb1, nb2, nb3; T, P and muB are calclated with EoS
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
  freezeoutSurfaceFile.open("output/surface.dat");

  //************************************************************************************\
  //* Fluid dynamics initialization
  //************************************************************************************/
    
  double t = t0;
  // Generate initial conditions
  setInitialConditions(latticeParams, initCondParams, hydroParams, rootDirectory);
  // Calculate conserved quantities
  setConservedVariables(t, latticeParams);
  // Slow modes for Hydro+
  setInitialConditionSlowModes(latticeParams, hydroParams);
  // Impose boundary conditions with ghost cells
  setGhostCells(q,e,p,u,latticeParams,rhob,alphaB,T,seq,eqPhiQ);

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
    
  FILE *fpan = fopen("output/AnalysisData.dat", "w");
    
  //************************************************************************************\
  //* loop over time steps
  //************************************************************************************/
    
  cout << "=============================================================="   << endl;
    
  for (int n = 1; n <= nt+1; ++n)
  {
    
    outputAnalysisa(n, t, fpan, latticeParams);
    //outputBaryonCP(t, outputDir, latticeParams);
      
    // copy variables back to host and write to disk
    if ((n-1) % FREQ == 0)
    {
      printf("n = %d:%d (t = %.3f),\t (e, p) = (%.3f, %.3f) [fm^-4],\t (rhob = %.3f ),\t (T = %.3f [GeV]),\n",
      n - 1, nt, t, e[sctr], p[sctr], rhob[sctr], effectiveTemperature(e[sctr], rhob[sctr]) * hbarc);
      outputDynamicalQuantities(t, outputDir, latticeParams);
      //outputHydroPlus(t, outputDir, latticeParams);
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
        callFOFinder(dim, start, nx, ny, nz, n, t0, dt, t, dx, dy, dz, lattice_spacing, freezeoutEnergyDensity, hyperCube4D, hyperCube3D,
                     energy_density_evoution, hydrodynamic_evoution, freezeoutSurfaceFile, fo_surf, FOFREQ);
    
    //if all cells are below freezeout temperature end hydro
    accumulator1 = 0;
    accumulator1 = checkForCellsAboveTc(nx, ny, nz, freezeoutEnergyDensity, e);
      
    if (accumulator1 == 0) accumulator2 += 1;
    if (accumulator2 >= FOFREQ+1 && (n > numberOfSourceFiles)) //only break once freezeout finder has had a chance to search/write to file
    {
      printf("\nAll cells have dropped below freezeout energy density\n");
      break;
    }
    
    //************************************************************************************\
    //* Evolution by Runge-Kutta
    //************************************************************************************/

    t1 = std::clock();

    //=======================================================
    //*  JET evolution
    //=======================================================
    if (JET)
    {
      //get the local fluid velocity and energy density/temperature and evolve jet momentum
      parton.energyLoss(nx, ny, nz, dt, dx, dy, dz, u->ut, u->ux, u->uy, u->un, e, rhob);
      //evolve the jet parton position
      parton.updatePosition(dt);
      //set hydro source terms
      setDynamicalSources(latticeParams, initCondParams, parton.dp_dtau, parton.position);
    }
    
    //=======================================================
    //*  dynamical sources
    //=======================================================
      
    // Read in source terms from particles
#ifdef DYNAMICAL_SOURCE
    if(initialConditionType==13){
          if(n <= numberOfSourceFiles) readInSource(n, latticeParams, initCondParams, hydroParams, rootDirectory);
          else if(n == numberOfSourceFiles + 1) zeroSource(latticeParams, initCondParams);
    }
#endif
      
    //=======================================================
    //*  RK-KT algorithm
    //=======================================================

    rungeKutta2(t, dt, q, Q, latticeParams, hydroParams);
    setCurrentConservedVariables();
      
    //=======================================================
    //*  Info
    //=======================================================
      
    t2 = std::clock();
    double delta_time = (t2 - t1) / (double)(CLOCKS_PER_SEC / 1000);
    if ((n-1) % FREQ == 0) printf("(Elapsed time: %.3f ms)\n",delta_time);
    totalTime+=delta_time;
    ++nsteps;

    t = t0 + n * dt;
  }
    
  printf("Average time/step: %.3f ms\n",totalTime/((double)nsteps));
    
  fclose(fpan);
  //************************************************************************************\
  //* Deallocate memory
  //************************************************************************************/
    
  freezeoutSurfaceFile.close();
    
  freeHostMemory();
  
  printf("freeHostMemory!\n");
    
  //Deallocate memory used for freezeout finding
  free4dArray(energy_density_evoution, FOFREQ+1, nx, ny);
  free5dArray(hydrodynamic_evoution, n_hydro_vars, FOFREQ+1, nx, ny);
  delete [] lattice_spacing;
    
  free4dArray(hyperCube4D, 2, 2, 2);
}
