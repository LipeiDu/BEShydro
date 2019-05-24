double ** calloc2dArray(double **array, int dim1, int dim2)
{
  array = (double **)calloc(dim1, sizeof(double *));
  for (int i1 = 0; i1 < dim1; i1++)
  {
    array[i1] = (double *)calloc(dim2, sizeof(double));
  }
  return array;
}

double *** calloc3dArray(double ***array, int dim1, int dim2, int dim3)
{
  array = (double ***)calloc(dim1, sizeof(double **));
  for (int i1 = 0; i1 < dim1; i1++)
  {
    array[i1] = (double **)calloc(dim2, sizeof(double *));
    for (int i2 = 0; i2 < dim2; i2++)
    {
      array[i1][i2] = (double *)calloc(dim3, sizeof(double));
    }
  }
  return array;
}

double **** calloc4dArray(double ****array, int dim1, int dim2, int dim3, int dim4)
{
  array = (double ****)calloc(dim1, sizeof(double ***));
  for (int i1 = 0; i1 < dim1; i1++)
  {
    array[i1] = (double ***)calloc(dim2, sizeof(double **));
    for (int i2 = 0; i2 < dim2; i2++)
    {
      array[i1][i2] = (double **)calloc(dim3, sizeof(double *));
      for (int i3 = 0; i3 < dim3; i3++)
      {
        array[i1][i2][i3] = (double *)calloc(dim4, sizeof(double));
      }
    }
  }
  return array;
}

double ***** calloc5dArray(double *****array, int dim1, int dim2, int dim3, int dim4, int dim5)
{
  array = (double *****)calloc(dim1, sizeof(double ****));
  for (int i1 = 0; i1 < dim1; i1++)
  {
    array[i1] = (double ****)calloc(dim2, sizeof(double ***));
    for (int i2 = 0; i2 < dim2; i2++)
    {
      array[i1][i2] = (double ***)calloc(dim3, sizeof(double **));
      for (int i3 = 0; i3 < dim3; i3++)
      {
        array[i1][i2][i3] = (double **)calloc(dim4, sizeof(double *));
        for (int i4 = 0; i4 < dim4; i4++)
        {
          array[i1][i2][i3][i4] = (double *)calloc(dim5, sizeof(double));
        }
      }
    }
  }
  return array;
}

void free2dArray(double **array, int dim1)
{
  for (int i1 = 0; i1 < dim1; i1++)
  {
    free(array[i1]);
  }
  free(array);
}

void free3dArray(double ***array, int dim1, int dim2)
{
  for (int i1 = 0; i1 < dim1; i1++)
  {
    for (int i2 = 0; i2 < dim2; i2++)
    {
      free(array[i1][i2]);
    }
    free(array[i1]);
  }
  free(array);
}
void free4dArray(double ****array, int dim1, int dim2, int dim3)
{
  for (int i1 = 0; i1 < dim1; i1++)
  {
    for (int i2 = 0; i2 < dim2; i2++)
    {
      for (int i3 = 0; i3 < dim3; i3++)
      {
        free(array[i1][i2][i3]);
      }
      free(array[i1][i2]);
    }
    free(array[i1]);
  }
  free(array);
}

void free5dArray(double *****array, int dim1, int dim2, int dim3, int dim4)
{
  for (int i1 = 0; i1 < dim1; i1++)
  {
    for (int i2 = 0; i2 < dim2; i2++)
    {
      for (int i3 = 0; i3 < dim3; i3++)
      {
        for (int i4 = 0; i4 < dim4; i4++)
        {
          free(array[i1][i2][i3][i4]);
        }
        free(array[i1][i2][i3]);
      }
      free(array[i1][i2]);
    }
    free(array[i1]);
  }
  free(array);
}
