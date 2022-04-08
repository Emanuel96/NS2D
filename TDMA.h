#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void TDMA_X( double **PC, double **R, double *A, double *B, double *C, double *D, int NPX, int NPY, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dt )
{
    int i, j, i1;

    int k;

    double a, b, c;

    //CAUTION HERE b and c negative

    //double A[NC], B[NC], C[NC], D[NC];

    for(k=1;k<=1;k++)
    {
        A[0] = 0;

        for(j=1;i<=NPY-2;j++)
        {
            C[0]=PC[0][j];

            for(i=1;i<=NPX-2;i++)
            {
                //Coefficients

                a = 2*(dt/(DX[i+1]*DX[i])+dt/(DY[j+1]*DY[j]));
                b = -1*dt/(DY[j+1]*DY[j]);
                c = -1*dt/(DX[i+1]*DX[i]);

                A[i]=-c;
                B[i]=-c;
                C[i]=-(b*PC[i][j+1]+b*PC[i][j-1]+R[i][j]);
                D[i]=a;

                //Recurrence Formula
                A[i]=A[i]/(D[i]-B[i]*A[i-1]);
                C[i]=(B[i]*C[i-1]+C[i])/(D[i]-B[i]*A[i-1]);
            }

            //New Phi's
            for(i1=NPX-2;i1>=1;i1--)
            {
                PC[i1][j]=A[i1]*PC[i1+1][j]+C[i1];
            }
        }
    }
}

void TDMA_Y( double **PC, double **R, double *A, double *B, double *C, double *D, int NPX, int NPY, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dt )
{
    int i, j, j1;

    int k;

    double a, b, c;

    //CAUTION HERE b and c negative

    //double A[NC], B[NC], C[NC], D[NC];

    for(k=1;k<=1;k++)
    {
        A[0] = 0;

        for(i=1;i<=NPX-2;i++)
        {
            C[0]=PC[i][0];

            for(j=1;j<=NPY-2;j++)
            {
                //Coefficients

                a = 2*(dt/(DX[i+1]*DX[i])+dt/(DY[j+1]*DY[j]));
                b = -1*dt/(DX[i+1]*DX[i]);
                c = -1*dt/(DY[j+1]*DY[j]);

                A[j]=-c;
                B[j]=-c;
                C[j]=-(b*PC[i+1][j]+b*PC[i-1][j]+R[i][j]);
                D[j]=a;

                //Recurrence Formula
                A[j]=A[j]/(D[j]-B[j]*A[j-1]);
                C[j]=(B[j]*C[j-1]+C[j])/(D[j]-B[j]*A[j-1]);
            }

            //New Phi's
            for(j1=NPY-2;j1>=1;j1--)
            {
                PC[i][j1]=A[j1]*PC[i][j1+1]+C[j1];
                //printf("\n%d %d : %f",i,j1,PC[i][j1]);
                //getchar();
            }
        }
    }
}
