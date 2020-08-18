#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h"
#include "../include/FreezeOut.h"
#include "cornelius-c++-1.3/cornelius.cpp"
#include <stdio.h>
#include <iostream>
#include <fstream>


//return a 4 dimensional linear interpolation inside the hypercube, given the values
//on the corners (a0000 through a1111) and edge lengths x0 through x3
double linearInterp4D(double x0, double x1, double x2, double x3, double a0000, double a1000, double a0100, double a0010, double a0001,
                      double a1100, double a1010, double a1001, double a0110, double a0101, double a0011,
                      double a1110, double a1101, double a0111, double a1011, double a1111)
{
  double result = 0;
  result = ((1-x0) * (1-x1) * (1-x2) * (1-x3) * a0000)
            + ((x0) * (1-x1) * (1-x2) * (1-x3) * a1000)
            + ((1-x0) * (x1) * (1-x2) * (1-x3) * a0100)
            + ((1-x0) * (1-x1) * (x2) * (1-x3) * a0010)
            + ((1-x0) * (1-x1) * (1-x2) * (x3) * a0001)
            + ((x0) * (x1) * (1-x2) * (1-x3) * a1100)
            + ((x0) * (1-x1) * (x2) * (1-x3) * a1010)
            + ((x0) * (1-x1) * (1-x2) * (x3) * a1001)
            + ((1-x0) * (x1) * (x2) * (1-x3) * a0110)
            + ((1-x0) * (x1) * (1-x2) * (x3) * a0101)
            + ((1-x0) * (1-x1) * (x2) * (x3) * a0011)
            + ((x0) * (x1) * (x2) * (1-x3) * a1110)
            + ((x0) * (x1) * (1-x2) * (x3) * a1101)
            + ((x0) * (1-x1) * (x2) * (x3) * a1011)
            + ((1-x0) * (x1) * (x2) * (x3) * a0111)
            + ((x0) * (x1) * (x2) * (x3) * a1111);

  return result;
}

double linearInterp3D(double x0, double x1, double x2, double a000, double a100, double a010, double a001, double a110, double a101, double a011, double a111)
{
  double result = 0;
  result = ((1-x0) * (1-x1) * (1-x2) * a000)
            + ((x0) * (1-x1) * (1-x2) * a100)
            + ((1-x0) * (x1) * (1-x2) * a010)
            + ((1-x0) * (1-x1) * (x2) * a001)
            + ((x0) * (x1) * (1-x2) * a110)
            + ((x0) * (1-x1) * (x2) * a101)
            + ((1-x0) * (x1) * (x2) * a011)
            + ((x0) * (x1) * (x2)  * a111);

  return result;
}

void swapAndSetHydroVariables(double ****energy_density_evoution, double *****hydrodynamic_evoution, CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int FOFREQ)
{
  for (int ix = 2; ix < nx+2; ix++)
  {
    for (int iy = 2; iy < ny+2; iy++)
    {
      for (int iz = 2; iz < nz+2; iz++)
      {
        int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4);
        //previous hydro variable values written to zeroth index
        energy_density_evoution[0][ix-2][iy-2][iz-2] = energy_density_evoution[FOFREQ][ix-2][iy-2][iz-2];
          
        //u0, u1, u2, u3, e, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, Pi, rhob, nb0, nb1, nb2, nb3;
        hydrodynamic_evoution[0][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[0][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[1][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[1][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[2][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[2][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[3][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[3][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[4][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[4][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[5][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[5][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[6][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[6][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[7][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[7][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[8][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[8][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[9][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[9][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[10][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[10][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[11][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[11][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[12][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[12][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[13][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[13][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[14][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[14][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[15][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[15][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[16][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[16][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[17][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[17][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[18][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[18][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[19][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[19][FOFREQ][ix-2][iy-2][iz-2];
        hydrodynamic_evoution[20][0][ix-2][iy-2][iz-2] = hydrodynamic_evoution[20][FOFREQ][ix-2][iy-2][iz-2];


        //current hydro variable values written to first index
        energy_density_evoution[1][ix-2][iy-2][iz-2] = (double)e[s];
        hydrodynamic_evoution[0][1][ix-2][iy-2][iz-2] = (double)(u->ut[s]);
        hydrodynamic_evoution[1][1][ix-2][iy-2][iz-2] = (double)(u->ux[s]);
        hydrodynamic_evoution[2][1][ix-2][iy-2][iz-2] = (double)(u->uy[s]);
        hydrodynamic_evoution[3][1][ix-2][iy-2][iz-2] = (double)(u->un[s]);
        hydrodynamic_evoution[4][1][ix-2][iy-2][iz-2] = (double)(e[s]);
#ifdef PIMUNU
        hydrodynamic_evoution[5][1][ix-2][iy-2][iz-2] = (double)(q->pitt[s]);
        hydrodynamic_evoution[6][1][ix-2][iy-2][iz-2] = (double)(q->pitx[s]);
        hydrodynamic_evoution[7][1][ix-2][iy-2][iz-2] = (double)(q->pity[s]);
        hydrodynamic_evoution[8][1][ix-2][iy-2][iz-2] = (double)(q->pitn[s]);
        hydrodynamic_evoution[9][1][ix-2][iy-2][iz-2] = (double)(q->pixx[s]);
        hydrodynamic_evoution[10][1][ix-2][iy-2][iz-2] = (double)(q->pixy[s]);
        hydrodynamic_evoution[11][1][ix-2][iy-2][iz-2] = (double)(q->pixn[s]);
        hydrodynamic_evoution[12][1][ix-2][iy-2][iz-2] = (double)(q->piyy[s]);
        hydrodynamic_evoution[13][1][ix-2][iy-2][iz-2] = (double)(q->piyn[s]);
        hydrodynamic_evoution[14][1][ix-2][iy-2][iz-2] = (double)(q->pinn[s]);
#endif
#ifdef PI
        hydrodynamic_evoution[15][1][ix-2][iy-2][iz-2] = (double)(q->Pi[s]);
#endif
#ifdef NBMU
        hydrodynamic_evoution[16][1][ix-2][iy-2][iz-2] = (double)(rhob[s]);
#endif
#ifdef VMU
        hydrodynamic_evoution[17][1][ix-2][iy-2][iz-2] = (double)(q->nbt[s]);
        hydrodynamic_evoution[18][1][ix-2][iy-2][iz-2] = (double)(q->nbx[s]);
        hydrodynamic_evoution[19][1][ix-2][iy-2][iz-2] = (double)(q->nby[s]);
        hydrodynamic_evoution[20][1][ix-2][iy-2][iz-2] = (double)(q->nbn[s]);
#endif
      }
    }
  }
}


void setHydroVariables(double ****energy_density_evoution, double *****hydrodynamic_evoution, CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int FOFREQ, int n)
{
  int nFO = n % FOFREQ;
  for (int ix = 2; ix < nx+2; ix++)
  {
    for (int iy = 2; iy < ny+2; iy++)
    {
      for (int iz = 2; iz < nz+2; iz++)
      {
        int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4);
        energy_density_evoution[nFO+1][ix-2][iy-2][iz-2] = (double)e[s];
        hydrodynamic_evoution[0][nFO+1][ix-2][iy-2][iz-2] = (double)(u->ut[s]);
        hydrodynamic_evoution[1][nFO+1][ix-2][iy-2][iz-2] = (double)(u->ux[s]);
        hydrodynamic_evoution[2][nFO+1][ix-2][iy-2][iz-2] = (double)(u->uy[s]);
        hydrodynamic_evoution[3][nFO+1][ix-2][iy-2][iz-2] = (double)(u->un[s]);
        hydrodynamic_evoution[4][nFO+1][ix-2][iy-2][iz-2] = (double)(e[s]);
#ifdef PIMUNU
        hydrodynamic_evoution[5][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pitt[s]);
        hydrodynamic_evoution[6][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pitx[s]);
        hydrodynamic_evoution[7][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pity[s]);
        hydrodynamic_evoution[8][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pitn[s]);
        hydrodynamic_evoution[9][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pixx[s]);
        hydrodynamic_evoution[10][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pixy[s]);
        hydrodynamic_evoution[11][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pixn[s]);
        hydrodynamic_evoution[12][nFO+1][ix-2][iy-2][iz-2] = (double)(q->piyy[s]);
        hydrodynamic_evoution[13][nFO+1][ix-2][iy-2][iz-2] = (double)(q->piyn[s]);
        hydrodynamic_evoution[14][nFO+1][ix-2][iy-2][iz-2] = (double)(q->pinn[s]);
#endif
#ifdef PI
        hydrodynamic_evoution[15][nFO+1][ix-2][iy-2][iz-2] = (double)(q->Pi[s]);
#endif
#ifdef NBMU
        hydrodynamic_evoution[16][nFO+1][ix-2][iy-2][iz-2] = (double)(rhob[s]);
#endif
#ifdef VMU
        hydrodynamic_evoution[17][nFO+1][ix-2][iy-2][iz-2] = (double)(q->nbt[s]);
        hydrodynamic_evoution[18][nFO+1][ix-2][iy-2][iz-2] = (double)(q->nbx[s]);
        hydrodynamic_evoution[19][nFO+1][ix-2][iy-2][iz-2] = (double)(q->nby[s]);
        hydrodynamic_evoution[20][nFO+1][ix-2][iy-2][iz-2] = (double)(q->nbn[s]);
#endif
      }
    }
  }
}


void writeEnergyDensityToHypercube4D(double ****hyperCube, double ****energy_density_evoution, int it, int ix, int iy, int iz)
{
  hyperCube[0][0][0][0] = energy_density_evoution[it][ix][iy][iz];
  hyperCube[1][0][0][0] = energy_density_evoution[it+1][ix][iy][iz];
  hyperCube[0][1][0][0] = energy_density_evoution[it][ix+1][iy][iz];
  hyperCube[0][0][1][0] = energy_density_evoution[it][ix][iy+1][iz];
  hyperCube[0][0][0][1] = energy_density_evoution[it][ix][iy][iz+1];
  hyperCube[1][1][0][0] = energy_density_evoution[it+1][ix+1][iy][iz];
  hyperCube[1][0][1][0] = energy_density_evoution[it+1][ix][iy+1][iz];
  hyperCube[1][0][0][1] = energy_density_evoution[it+1][ix][iy][iz+1];
  hyperCube[0][1][1][0] = energy_density_evoution[it][ix+1][iy+1][iz];
  hyperCube[0][1][0][1] = energy_density_evoution[it][ix+1][iy][iz+1];
  hyperCube[0][0][1][1] = energy_density_evoution[it][ix][iy+1][iz+1];
  hyperCube[1][1][1][0] = energy_density_evoution[it+1][ix+1][iy+1][iz];
  hyperCube[1][1][0][1] = energy_density_evoution[it+1][ix+1][iy][iz+1];
  hyperCube[1][0][1][1] = energy_density_evoution[it+1][ix][iy+1][iz+1];
  hyperCube[0][1][1][1] = energy_density_evoution[it][ix+1][iy+1][iz+1];
  hyperCube[1][1][1][1] = energy_density_evoution[it+1][ix+1][iy+1][iz+1];
}

void writeEnergyDensityToHypercube3D(double ***hyperCube, double ****energy_density_evoution, int it, int ix, int iy)
{
  hyperCube[0][0][0] = energy_density_evoution[it][ix][iy][0];
  hyperCube[1][0][0] = energy_density_evoution[it+1][ix][iy][0];
  hyperCube[0][1][0] = energy_density_evoution[it][ix+1][iy][0];
  hyperCube[0][0][1] = energy_density_evoution[it][ix][iy+1][0];
  hyperCube[1][1][0] = energy_density_evoution[it+1][ix+1][iy][0];
  hyperCube[1][0][1] = energy_density_evoution[it+1][ix][iy+1][0];
  hyperCube[0][1][1] = energy_density_evoution[it][ix+1][iy+1][0];
  hyperCube[1][1][1] = energy_density_evoution[it+1][ix+1][iy+1][0];
}

double interpolateVariable4D(double *****hydrodynamic_evoution, int ivar, int it, int ix, int iy, int iz, double tau_frac, double x_frac, double y_frac, double z_frac)
{
  double result = linearInterp4D(tau_frac, x_frac, y_frac, z_frac,
                                 hydrodynamic_evoution[ivar][it][ix][iy][iz], hydrodynamic_evoution[ivar][it+1][ix][iy][iz],
                                 hydrodynamic_evoution[ivar][it][ix+1][iy][iz], hydrodynamic_evoution[ivar][it][ix][iy+1][iz],
                                 hydrodynamic_evoution[ivar][it][ix][iy][iz+1],hydrodynamic_evoution[ivar][it+1][ix+1][iy][iz],
                                 hydrodynamic_evoution[ivar][it+1][ix][iy+1][iz],hydrodynamic_evoution[ivar][it+1][ix][iy][iz+1],
                                 hydrodynamic_evoution[ivar][it][ix+1][iy+1][iz], hydrodynamic_evoution[ivar][it][ix+1][iy][iz+1],
                                 hydrodynamic_evoution[ivar][it][ix][iy+1][iz+1], hydrodynamic_evoution[ivar][it+1][ix+1][iy+1][iz],
                                 hydrodynamic_evoution[ivar][it+1][ix+1][iy][iz+1], hydrodynamic_evoution[ivar][it][ix+1][iy+1][iz+1],
                                 hydrodynamic_evoution[ivar][it+1][ix][iy+1][iz+1], hydrodynamic_evoution[ivar][it+1][ix+1][iy+1][iz+1]);
    return result;
}

double interpolateVariable3D(double *****hydrodynamic_evoution, int ivar, int it, int ix, int iy, double tau_frac, double x_frac, double y_frac)
{
  double result = linearInterp3D(tau_frac, x_frac, y_frac,
                                 hydrodynamic_evoution[ivar][it][ix][iy][0], hydrodynamic_evoution[ivar][it+1][ix][iy][0],
                                 hydrodynamic_evoution[ivar][it][ix+1][iy][0], hydrodynamic_evoution[ivar][it][ix][iy+1][0],
                                 hydrodynamic_evoution[ivar][it+1][ix+1][iy][0], hydrodynamic_evoution[ivar][it+1][ix][iy+1][0],
                                 hydrodynamic_evoution[ivar][it][ix+1][iy+1][0], hydrodynamic_evoution[ivar][it+1][ix+1][iy+1][0]);
    return result;
}

void callFOFinder(int dim, int start, int nx, int ny, int nz, int n, double t0, double dt, double t, double dx, double dy, double dz, double *lattice_spacing, double freezeoutEnergyDensity, double ****hyperCube4D, double ***hyperCube3D, double ****energy_density_evoution, double *****hydrodynamic_evoution, std::ofstream& freezeoutSurfaceFile, std::vector<FO_Element>& fo_surf, int FOFREQ)
{
    
  Cornelius cor;
  cor.init(dim, freezeoutEnergyDensity, lattice_spacing);
    
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
            double cell_tau;

            if (n <= FOFREQ)
                cell_tau = t0 + ( (double)(n - FOFREQ + (it-1) ) )* dt; //check if this is the correct time!
            else
                cell_tau = t0 + ((double)(n - FOFREQ + it)) * dt; //check if this is the correct time!
              
            double cell_x = (double)ix * dx  - (((double)(nx-1)) / 2.0 * dx);
            double cell_y = (double)iy * dy  - (((double)(ny-1)) / 2.0 * dy);
            double cell_z = (double)iz * dz  - (((double)(nz-1)) / 2.0 * dz);

            double tau_frac = cor.get_centroid_elem(i,0) / lattice_spacing[0];
            double x_frac = cor.get_centroid_elem(i,1) / lattice_spacing[1];
            double y_frac = cor.get_centroid_elem(i,2) / lattice_spacing[2];
            double z_frac;
            if (dim == 4) z_frac = cor.get_centroid_elem(i,3) / lattice_spacing[3];
            else z_frac = 0.0;

            //first write the contravariant position vector
            freezeoutSurfaceFile << cor.get_centroid_elem(i,0) + cell_tau << "\t ";
            freezeoutSurfaceFile << cor.get_centroid_elem(i,1) + cell_x << "\t ";
            freezeoutSurfaceFile << cor.get_centroid_elem(i,2) + cell_y << "\t ";
            if (dim == 4) freezeoutSurfaceFile << cor.get_centroid_elem(i,3) + cell_z << "\t ";
            else freezeoutSurfaceFile << cell_z << "\t ";
            //then the contravariant surface normal element; note jacobian factors of tau for milne coordinates
            freezeoutSurfaceFile << t * cor.get_normal_elem(i,0) << "\t ";
            freezeoutSurfaceFile << t * cor.get_normal_elem(i,1) << "\t ";
            freezeoutSurfaceFile << t * cor.get_normal_elem(i,2) << "\t ";
            if (dim == 4) freezeoutSurfaceFile << t * cor.get_normal_elem(i,3) << "\t ";
            else freezeoutSurfaceFile << 0.0 << "\t ";
            //write all the necessary hydro dynamic variables by first performing linear interpolation from values at
            //corners of hypercube
              
              
            //************************************************************************************\
            //* for 3+1D
            //************************************************************************************/
            if (dim == 4)
            {
              //first write the contravariant flow velocity
              for (int ivar = 0; ivar < dim; ivar++)
              {
                temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                freezeoutSurfaceFile << temp << "\t ";
              }
              //write the energy density
              temp = interpolateVariable4D(hydrodynamic_evoution, 4, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
              freezeoutSurfaceFile << temp << "\t "; //note : iSpectra reads in file in fm^x units e.g. energy density should be written in fm^-4
              //baryon density
              double temp2 = interpolateVariable4D(hydrodynamic_evoution, 16, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
              //the temperature !this needs to be checked
              double temp3 = effectiveTemperature(temp,temp2);
              freezeoutSurfaceFile << temp3 << "\t ";
              //the baryon chemical potential
              freezeoutSurfaceFile << chemicalPotentialOverT(temp,temp2)*temp3 << "\t ";
              //the thermal pressure
              freezeoutSurfaceFile << equilibriumPressure(temp,temp2) << "\t ";
              //write ten components of pi_(mu,nu) shear viscous tensor
              for (int ivar = 5; ivar < 15; ivar++)
              {
                temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                freezeoutSurfaceFile << temp << "\t ";
              }
              //write the bulk pressure Pi
              temp = interpolateVariable4D(hydrodynamic_evoution, 15, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
              freezeoutSurfaceFile << temp << "\t ";
              //write baryon density
              freezeoutSurfaceFile << temp2 << "\t ";
              //write the baryon diffusion
              for (int ivar = 17; ivar < 21; ivar++)
              {
                temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                freezeoutSurfaceFile << temp << "\t ";
              }
              //start a new line
              freezeoutSurfaceFile << endl;
            }
            //************************************************************************************\
            //* for 2+1D
            //************************************************************************************/
            else
            {
                //first write the contravariant flow velocity
                for (int ivar = 0; ivar < 4; ivar++)
                {
                  temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                  freezeoutSurfaceFile << setprecision(5) << setw(15) << temp ;//<< "\t ";
                }
                //write the energy density
                temp = interpolateVariable3D(hydrodynamic_evoution, 4, it, ix, iy, tau_frac, x_frac, y_frac);
                freezeoutSurfaceFile << setprecision(5) << setw(15) << temp ;//<< "\t "; //note units of fm^-4 appropriate for iSpectra reading
                //baryon density
                double temp2 = interpolateVariable4D(hydrodynamic_evoution, 16, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                //the temperature !this needs to be checked
                double temp3 = effectiveTemperature(temp,temp2);
                freezeoutSurfaceFile << setprecision(5) << setw(15) << temp3 ;//<< "\t ";
                //calculate chemical potential
                double temp4 = chemicalPotentialOverT(temp,temp2)*temp3;
                //the thermal pressure
                freezeoutSurfaceFile << setprecision(5) << setw(15) << equilibriumPressure(temp,temp2) ;//<< "\t ";
                //write ten components of pi_(mu,nu) shear viscous tensor
                for (int ivar = 5; ivar < 15; ivar++)
                {
                  temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                  freezeoutSurfaceFile << setprecision(5) << setw(15) << temp ;//<< "\t ";
                }
                //write the bulk pressure Pi
                temp = interpolateVariable3D(hydrodynamic_evoution, 15, it, ix, iy, tau_frac, x_frac, y_frac);
                freezeoutSurfaceFile << setprecision(5) << setw(15) << temp ;//<< "\t ";
                //write the baryon chemical potential
                freezeoutSurfaceFile << setprecision(5) << setw(15) << temp4 ;//<< "\t ";
                //write baryon density
                freezeoutSurfaceFile << setprecision(5) << setw(15) << temp2 ;//<< "\t ";
                //write the baryon diffusion
                for (int ivar = 17; ivar < 21; ivar++)
                {
                  temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                  freezeoutSurfaceFile << setprecision(5) << setw(15) << temp ;//<< "\t ";
                }
                //start a new line
                freezeoutSurfaceFile << setprecision(5) << setw(15) << endl;
            }
          }
        }
      }
    }
  }
}

//returns the number of cells with T > T_c
int checkForCellsAboveTc(int nx, int ny, int nz, double freezeoutEnergyDensity, PRECISION *e)
{
 int accumulator = 0;

 for (int ix = 2; ix < nx+2; ix++)
  {
   for (int iy = 2; iy < ny+2; iy++)
    {
     for (int iz = 2; iz < nz+2; iz++)
      {
        int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4);
        if (e[s] > freezeoutEnergyDensity) accumulator += 1;
      }
    }
  }

 return accumulator;
}

