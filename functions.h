#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double Maximum_Element_Matrix_NO_BOUNDARIES(double **X, int lin ,int col )
{
    int i,j;
    double maior = X[0][0];

    for ( i = 0 ; i < lin ; i++)
    {
        for ( j = 0 ; j < col ; j++)
        {
            if ( X[i][j] >= maior )
            {
                maior = X[i][j];
            }
        }

    }
    return maior;
}

void Transform_Absolute_Element_Matrix(double **X ,int lin ,int col )
{
    int i,j;

    for ( i = 0 ; i < lin ; i++)
    {
        for ( j = 0 ; j < col ; j++)
        {
            X[i][j]=fabs(X[i][j]);
        }
    }
}

