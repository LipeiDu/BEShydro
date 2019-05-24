/**
 *
 * Simple program which demostrates how one should use cornelius++
 * from cornelius.cpp. Makes squares, cubes and hypercubes
 * with random corner values from 0-1 and finds the surface elements
 * for the value 0.5 and writes them into files.
 *
 * Last update 24.07.2012 Hannu Holopainen
 *
 */

#include <iostream>
#include <fstream>
#include "cornelius.cpp"

using namespace std;

void example_2d()
{
  ofstream output;
  output.open("output_2d.dat");
  //Let's initialize cornelius first
  int dim = 2;
  double *dx = new double[dim];
  dx[0] = 1;
  dx[1] = 1;
  Cornelius cor;
  cor.init(dim,0.5,dx);
  //Allocating memory for corner points
  double **cube = new double*[2];
  for (int i1=0; i1 < 2; i1++) {
    cube[i1] = new double[2];
  }
  //Let's make 10 000 random cubes as a test case
  for (int iter=0; iter < 10000; iter++) {
    for (int i1=0; i1 < 2; i1++) {
      for (int i2=0; i2 < 2; i2++) {
        cube[i1][i2] = rand()/double(RAND_MAX);
      }
    }
    cor.find_surface_2d(cube);
    for (int i=0; i < cor.get_Nelements(); i++) {
      //First we write centroid
      for (int j=0; j < dim; j++) {
        output << cor.get_centroid_elem(i,j) << " ";
      }
      //Then normal
      for (int j=0; j < dim; j++) {
        output << cor.get_normal_elem(i,j) << " ";
      }
      output << endl;
    }
  }
  //All done, let's free memory
  for (int i1=0; i1 < 2; i1++) {
    delete[] cube[i1];
  }
  delete[] cube;
  output.close();
}

void example_3d()
{
  ofstream output;
  output.open("output_3d.dat");
  //Let's initialize cornelius first
  int dim = 3;
  double *dx = new double[dim];
  dx[0] = 1;
  dx[1] = 1;
  dx[2] = 1;
  Cornelius cor;
  cor.init(dim,0.5,dx);
  //Allocating memory for corner points
  double ***cube = new double**[2];
  for (int i1=0; i1 < 2; i1++) {
    cube[i1] = new double*[2];
    for (int i2=0; i2 < 2; i2++) {
      cube[i1][i2] = new double[2];
    }
  }
  //Let's make 10 000 random cubes as a test case
  for (int iter=0; iter < 10000; iter++) {
    for (int i1=0; i1 < 2; i1++) {
      for (int i2=0; i2 < 2; i2++) {
        for (int i3=0; i3 < 2; i3++) {
          cube[i1][i2][i3] = rand()/double(RAND_MAX);
        }
      }
    }
    cor.find_surface_3d(cube);
    for (int i=0; i < cor.get_Nelements(); i++) {
      //First we write centroid
      for (int j=0; j < dim; j++) {
        output << cor.get_centroid_elem(i,j) << " ";
      }
      //Then normal
      for (int j=0; j < dim; j++) {
        output << cor.get_normal_elem(i,j) << " ";
      }
      output << endl;
    }
  }
  //All done, let's free memory
  for (int i1=0; i1 < 2; i1++) {
    for (int i2=0; i2 < 2; i2++) {
      delete[] cube[i1][i2];
    }
    delete[] cube[i1];
  }
  delete[] cube;
  output.close();
}

void example_4d()
{
  ofstream output;
  output.open("output_4d.dat");
  //Let's initialize cornelius first
  int dim = 4;
  double *dx = new double[dim];
  dx[0] = 1;
  dx[1] = 1;
  dx[2] = 1;
  dx[3] = 1;
  Cornelius cor;
  cor.init(dim,0.5,dx);
  //Allocating memory for corner points
  double ****cube = new double***[2];
  for (int i1=0; i1 < 2; i1++) {
    cube[i1] = new double**[2];
    for (int i2=0; i2 < 2; i2++) {
      cube[i1][i2] = new double*[2];
      for (int i3=0; i3 < 2; i3++) {
        cube[i1][i2][i3] = new double[2];
      }
    }
  }
  //Let's make 10 000 random cubes as a test case
  for (int iter=0; iter < 10000; iter++) {
    for (int i1=0; i1 < 2; i1++) {
      for (int i2=0; i2 < 2; i2++) {
        for (int i3=0; i3 < 2; i3++) {
          for (int i4=0; i4 < 2; i4++) {
            cube[i1][i2][i3][i4] = rand()/double(RAND_MAX);
          }
        }
      }
    }
    cor.find_surface_4d(cube);
    for (int i=0; i < cor.get_Nelements(); i++) {
      //First we write centroid
      for (int j=0; j < dim; j++) {
        output << cor.get_centroid_elem(i,j) << " ";
      }
      //Then normal
      for (int j=0; j < dim; j++) {
        output << cor.get_normal_elem(i,j) << " ";
      }
      output << endl;
    }
  }
  //All done, let's free memory
  for (int i1=0; i1 < 2; i1++) {
    for (int i2=0; i2 < 2; i2++) {
      for (int i3=0; i3 < 2; i3++) {
        delete[] cube[i1][i2][i3];
      }
      delete[] cube[i1][i2];
    }
    delete[] cube[i1];
  }
  delete[] cube;
  output.close();
}

int main()
{
  example_2d(); 
  example_3d(); 
  example_4d(); 
  return 1;
}

