/*
 * FileIO.c
 *
 *  Created on: Oct 24, 2015
 *      Author: bazow
 */
#include <stdlib.h>
#include <stdio.h>

#include "edu/osu/rhic/harness/io/FileIO.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

void output(const PRECISION * const var, double t, const char *pathToOutDir, const char *name, void * latticeParams) {
	FILE *fp;
	char fname[255];
	sprintf(fname, "%s/%s_%.3f.dat", pathToOutDir, name, t);
	fp=fopen(fname, "w");

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
/*
	for(i = 2; i < nx+2; ++i) {
		x = (i-2 - (nx-1)/2.)*dx;
		for(j = 2; j < ny+2; ++j) {
			y = (j-2 - (ny-1)/2.)*dy;
			offset = (nz+4) * (j + (ny+4) * i);
			for(k = 2; k < nz+2; ++k) {
				s = k + offset;
				z = (k-2 - (nz-1)/2.)*dz;
				fprintf(fp, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,var[s]);
			}
		}
	}
*/
	for(k = 2; k < nz+2; ++k) {
		z = (k-2 - (nz-1)/2.)*dz;
		for(j = 2; j < ny+2; ++j) {
			y = (j-2 - (ny-1)/2.)*dy;
			for(i = 2; i < nx+2; ++i) {
				x = (i-2 - (nx-1)/2.)*dx;
				s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                //if(x==0&&y==0)
				fprintf(fp, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,var[s]);
			}
		}
	}

	fclose(fp);
}
