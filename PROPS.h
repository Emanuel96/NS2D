#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void PROPS( double **T, double **P, double vis, double **VIS, double **VIS_T, double **RHO, double **CP, double **K, int NPX, int NPY, double **U, double **V, double *X, double *Y, double *DX, double *DY, double *DXU, double *DYU, double *DXV, double *DYV, double **I_x, double **I_y, double **Re, double dts )
{
    int i, j;

    //double Re_c = 20;

    for(i=0;i<=NPX-1;i++)
    {
        for(j=0;j<=NPY-1;j++)
        {
            //TURBULENCE
            //I_x[i][j] = 0.05/(1+exp(20*(Re[i][j]-Re_c)));
            //I_y[i][j] = 0.05/(1+exp(20*(Re[i][j]-Re_c)));

            if ( i == 0 || i == NPX || j == 0 || j == NPY )
                VIS_T[i][j] = 0;
            else
                VIS_T[i][j] = 1e-4*RHO[i][j]*fabs((V[i+1][j]-V[i][j])/DXV[i]-(U[i][j+1]-U[i][j])/DYU[j]);

            //DENSITY
            RHO[i][j] = P[i][j]/(287*T[i][j]);//RHO[i][j];

            //VISCOSITY
            VIS[i][j] = vis*pow(T[i][j]/273.11,3/2)*(273.11+110.56)/(T[i][j]+110.56);
        }
    }
}
