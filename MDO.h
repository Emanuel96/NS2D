#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double error_old = 0;
double error;
double integral = 0;
double derivative;

double PID( double **U, double **V, int NPX, int NPY, double dt, double U_in )
{
    double SP = 2*U_in;
    double K_P = 0.01;
    double T_I = 5.0;
    double T_D = 2.5;

    double MV;

    int j;

    double sum = 0;
    double n = 0;

    for( j=0; j<NPY; j++ )
    {
        if ( U[NPX-1][j] > 0 )
        {
            sum = sum + U[NPX-1][j];
            n = n + 1;
        }
    }

    MV = sum/n;

    error = SP - MV;
    integral = integral + error*dt;
    derivative = (error-error_old)/dt;
    error_old=error;

    //printf("\nVEL: %lf\n",K_P*error+K_I*integral+K_D*derivative);

    return K_P*(error+1/T_I*integral+1/T_D*derivative);
}
