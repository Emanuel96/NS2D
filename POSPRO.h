#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void LAW_WALL(double vis, double rho, int NPX, int NPY, double **U, double **V, double *X, double *Y, double *DX, double *DY, double *DXU, double *DYU, double *DXV, double *DYV )
{
    int i, j;

    char name[50];

    double u_tau;
    double y_plus, u_plus;
    double HEIGHT_BOUNDARY_LAYER = 0.4;

    for(i=0; i<=NPX-1; i++)
    {
        if ( i == 50 || i == 60 || i == 70 || i == 80 || i == 90 )
        {
            FILE* file_law;
            sprintf(name,"./Results/%d.dat",i);
            file_law = fopen(name,"w");
            fprintf(file_law,"VARIABLES=y<sup>+</sup>,u<sup>+</sup>\n");

            u_tau = sqrt(vis*(U[i][1]/DYU[1])/rho);

            for(j=0; j<NPY; j++)
            {
                if ( Y[j] <= HEIGHT_BOUNDARY_LAYER )
                {
                    y_plus = Y[j]*u_tau*rho/vis;
                    u_plus = U[i][j]/u_tau;
                    fprintf(file_law,"%lf %lf\n",y_plus,u_plus);
                }
            }
            fclose(file_law);
        }
    }
}

//    FILE* file_C;
//    file_C = fopen("./Results/control.dat","w");
//
//    //Control point
//    //fprintf(file_C,"%lf %lf %lf\n",t*dts,U[40][25],V[40][25]);
//
//    fclose(file_C);
