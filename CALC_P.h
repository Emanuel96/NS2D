#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void CALC_P( double **P, double **PC, int NPX, int NPY, double urf )
{
    int i, j;

    for(i=1;i<=NPX-1;i++)
    {
        for(j=1;j<=NPY-1;j++)
        {
            P[i][j]=P[i][j]+urf*PC[i][j];
            PC[i][j]=0;
        }
    }
}
