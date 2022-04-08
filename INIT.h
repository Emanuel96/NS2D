#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void INIT( double **X, int NPX, int NPY, double value )
{
    int i, j;

    for(i=0;i<NPX;i++)
    {
        for(j=0;j<NPY;j++)
        {
            X[i][j] = value;
        }
    }
}
