#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define gx 0
#define gy 0

double random()
{
    return (double)rand()/RAND_MAX;
}

void CONTINUITY( double **U, double **V, double **RHO, double **R, int NPX, int NPY, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dt)
{
    int i, j;

    for(i=0;i<=NPX-1;i++)
    {
        for(j=0;j<=NPY-1;j++)
        {
            R[i][j] = (RHO[i][j]*(U[i+1][j]-U[i][j]))/DXU[i+1]+(RHO[i][j]*(V[i][j+1]-V[i][j]))/DYV[j+1];
        }
    }
}

void CALC_U_C( double **U, double **V, double **PC, double **RHO, int NPX, int NPY, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dt )
{
    int i, j;
    double p;
    double rho;

    for(i=2;i<=NPX-2;i++)
    {
        for(j=1;j<=NPY-2;j++)
        {
            rho = 0.5*(RHO[i][j]+RHO[i-1][j]);
            p = -1/DX[i]*(PC[i][j]-PC[i-1][j]);
            U[i][j] = U[i][j] + 1/rho*dt*p;
        }
    }
}

void CALC_V_C( double **U, double **V, double **PC, double **RHO, int NPX, int NPY, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dt )
{
    int i, j;
    double p;
    double rho;

    for(i=1;i<=NPX-2;i++)
    {
        for(j=2;j<=NPY-2;j++)
        {
            rho = 0.5*(RHO[i][j]+RHO[i][j-1]);
            p = -1/DY[j]*(PC[i][j]-PC[i][j-1]);
            V[i][j] = V[i][j] + 1/rho*dt*p;
        }
    }
}

void CALC_U( double **U, double **V, double **P, int NPX, int NPY, double *Y, double *YV, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dt, double **RHO, double **VIS, double **VIS_T )
{
    int i, j;
    double c1, c2, d, p, g, t;
    double v1, v2, v3, v4, vd, vu;
    double rho1, rho2, rho_u, rho_d, rho;

    for(i=2;i<=NPX-2;i++)
    {
        for(j=1;j<=NPY-2;j++)
        {
            rho_u = 0.5*(RHO[i+1][j]+RHO[i][j]);
            rho_d = 0.5*(RHO[i-1][j]+RHO[i-2][j]);

            c1 = -(U[i+1][j]*U[i+1][j]*rho_u-U[i-1][j]*U[i-1][j]*rho_d)/(DXU[i+1]+DXU[i]);

            v1 = 0.5*(V[i-1][j]+V[i][j]);
            v2 = 0.5*(V[i-1][j-1]+V[i][j-1]);
            v3 = 0.5*(V[i][j+2]+V[i-1][j+2]);
            v4 = 0.5*(V[i][j+1]+V[i-1][j+1]);

            rho1 = 0.5*(RHO[i-1][j+1]+RHO[i][j+1]);
            rho2 = 0.5*(RHO[i-1][j-1]+RHO[i][j-1]);

            vu = ((Y[j+1]-YV[j+1])/(YV[j+2]-YV[j+1]))*v3+(1-(Y[j+1]-YV[j+1])/(YV[j+2]-YV[j+1]))*v4;
            vd = ((Y[j-1]-YV[j-1])/(YV[j]-YV[j-1]))*v1+(1-(Y[j-1]-YV[j-1])/(YV[j]-YV[j-1]))*v2;

            c2 = -(U[i][j+1]*vu*rho1-U[i][j-1]*vd*rho2)/(DYU[j+1]+DYU[j]);

            p = -1/DX[i]*(P[i][j]-P[i-1][j]);

            d = VIS[i][j]*((U[i+1][j]-2*U[i][j]+U[i-1][j])/(DXU[i+1]*DXU[i])+(U[i][j+1]-2*U[i][j]+U[i][j-1])/(DYU[j+1]*DYU[j]));

            t = VIS_T[i][j]*((U[i+1][j]-2*U[i][j]+U[i-1][j])/(DXU[i+1]*DXU[i])+(U[i][j+1]-2*U[i][j]+U[i][j-1])/(DYU[j+1]*DYU[j]));

            g = gx;

            rho = 0.5*(RHO[i-1][j]+RHO[i][j]);

            U[i][j] = U[i][j] + 1/rho*dt*(c1 + c2 + p + d + g + t);
        }
    }
}

void CALC_V( double **U, double **V, double **P, int NPX, int NPY, double *X, double *XU, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dt, double **RHO, double **VIS, double **VIS_T )
{
    int i, j;
    double c1, c2, d, p, g, t;
    double u1, u2, u3, u4, ul, ur;
    double rho1, rho2, rho_u, rho_d, rho;

    for(i=1;i<=NPX-2;i++)
    {
        for(j=2;j<=NPY-2;j++)
        {
            rho_u = 0.5*(RHO[i][j+1]+RHO[i][j]);
            rho_d = 0.5*(RHO[i][j-1]+RHO[i][j-2]);

            c1 = -(V[i][j+1]*V[i][j+1]*rho_u-V[i][j-1]*V[i][j-1]*rho_d)/(DYV[j+1]+DYV[j]);

            u1 = 0.5*(U[i-1][j-1]+U[i-1][j]);
            u2 = 0.5*(U[i][j-1]+U[i][j]);
            u3 = 0.5*(U[i+1][j]+U[i+1][j-1]);
            u4 = 0.5*(U[i+2][j]+U[i+2][j-1]);

            rho1 = 0.5*(RHO[i+1][j]+RHO[i+1][j-1]);
            rho2 = 0.5*(RHO[i-1][j-1]+RHO[i-1][j]);

            ul = (X[i-1]-XU[i-1])/(XU[i]-XU[i-1])*u2+(1-(X[i-1]-XU[i-1])/(XU[i]-XU[i-1]))*u1;
            ur = (X[i+1]-XU[i+1])/(XU[i+2]-XU[i+1])*u4+(1-(X[i+1]-XU[i+1])/(XU[i+2]-XU[i+1]))*u3;

            c2 = -(V[i+1][j]*ur*rho1-V[i-1][j]*ul*rho2)/(DXV[i+1]+DXV[i]);

            p = -1/DY[j]*(P[i][j]-P[i][j-1]);

            d = VIS[i][j]*((V[i][j+1]-2*V[i][j]+V[i][j-1])/(DYV[j+1]*DYV[j])+(V[i+1][j]-2*V[i][j]+V[i-1][j])/(DXV[i+1]*DXV[i]));

            t = VIS_T[i][j]*((V[i][j+1]-2*V[i][j]+V[i][j-1])/(DYV[j+1]*DYV[j])+(V[i+1][j]-2*V[i][j]+V[i-1][j])/(DXV[i+1]*DXV[i]));

            g = gy;

            rho = 0.5*(RHO[i][j-1]+RHO[i][j]);

            V[i][j] = V[i][j] + 1/rho*dt*(c1 + c2 + p + d + g + t);
        }
    }
}

void CALC_U_t( double **U_NEW, double **P, double **U, double **V, int NPX, int NPY, double *Y, double *YV, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dts, double **I_x, double **RHO, double **VIS, double **VIS_T )
{
    int i, j;
    double c1, c2, d, p, g, t;
    double v1, v2, v3, v4, vd, vu;
    double rho1, rho2, rho_u, rho_d, rho;

    for(i=2;i<=NPX-2;i++)
    {
        for(j=1;j<=NPY-2;j++)
        {
            rho_u = 0.5*(RHO[i+1][j]+RHO[i][j]);
            rho_d = 0.5*(RHO[i-1][j]+RHO[i-2][j]);

            c1 = -(U[i+1][j]*U[i+1][j]*rho_u-U[i-1][j]*U[i-1][j]*rho_d)/(DXU[i+1]+DXU[i]);

            v1 = 0.5*(V[i-1][j]+V[i][j]);
            v2 = 0.5*(V[i-1][j-1]+V[i][j-1]);
            v3 = 0.5*(V[i][j+2]+V[i-1][j+2]);
            v4 = 0.5*(V[i][j+1]+V[i-1][j+1]);

            rho1 = 0.5*(RHO[i-1][j+1]+RHO[i][j+1]);
            rho2 = 0.5*(RHO[i-1][j-1]+RHO[i][j-1]);

            vu = ((Y[j+1]-YV[j+1])/(YV[j+2]-YV[j+1]))*v3+(1-(Y[j+1]-YV[j+1])/(YV[j+2]-YV[j+1]))*v4;
            vd = ((Y[j-1]-YV[j-1])/(YV[j]-YV[j-1]))*v1+(1-(Y[j-1]-YV[j-1])/(YV[j]-YV[j-1]))*v2;

            c2 = -(U[i][j+1]*vu*rho1-U[i][j-1]*vd*rho2)/(DYU[j+1]+DYU[j]);

            p = -1/DX[i]*(P[i][j]-P[i-1][j]);

            d = VIS[i][j]*((U[i+1][j]-2*U[i][j]+U[i-1][j])/(DXU[i+1]*DXU[i])+(U[i][j+1]-2*U[i][j]+U[i][j-1])/(DYU[j+1]*DYU[j]));

            t = VIS_T[i][j]*((U[i+1][j]-2*U[i][j]+U[i-1][j])/(DXU[i+1]*DXU[i])+(U[i][j+1]-2*U[i][j]+U[i][j-1])/(DYU[j+1]*DYU[j]));

            g = gx;

            rho = 0.5*(RHO[i-1][j]+RHO[i][j]);

            U_NEW[i][j] = U[i][j] + 1/rho*dts*(c1 + c2 + p + d + g + t) + U[i][j]*(2*random()-1)*I_x[i][j];
            U[i][j] = U_NEW[i][j];
        }
    }
}

void CALC_V_t( double **V_NEW, double **P, double **U, double **V, int NPX, int NPY, double *X, double *XU, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dts, double **I_y, double **RHO, double **VIS, double **VIS_T )
{
    int i, j;
    double c1, c2, d, p, g, t;
    double u1, u2, u3, u4, ul, ur;
    double rho1, rho2, rho_u, rho_d, rho;

    for(i=1;i<=NPX-2;i++)
    {
        for(j=2;j<=NPY-2;j++)
        {
            rho_u = 0.5*(RHO[i][j+1]+RHO[i][j]);
            rho_d = 0.5*(RHO[i][j-1]+RHO[i][j-2]);

            c1 = -(V[i][j+1]*V[i][j+1]*rho_u-V[i][j-1]*V[i][j-1]*rho_d)/(DYV[j+1]+DYV[j]);

            u1 = 0.5*(U[i-1][j-1]+U[i-1][j]);
            u2 = 0.5*(U[i][j-1]+U[i][j]);
            u3 = 0.5*(U[i+1][j]+U[i+1][j-1]);
            u4 = 0.5*(U[i+2][j]+U[i+2][j-1]);

            rho1 = 0.5*(RHO[i+1][j]+RHO[i+1][j-1]);
            rho2 = 0.5*(RHO[i-1][j-1]+RHO[i-1][j]);

            ul = (X[i-1]-XU[i-1])/(XU[i]-XU[i-1])*u2+(1-(X[i-1]-XU[i-1])/(XU[i]-XU[i-1]))*u1;
            ur = (X[i+1]-XU[i+1])/(XU[i+2]-XU[i+1])*u4+(1-(X[i+1]-XU[i+1])/(XU[i+2]-XU[i+1]))*u3;

            c2 = -(V[i+1][j]*ur*rho1-V[i-1][j]*ul*rho2)/(DXV[i+1]+DXV[i]);

            p = -1/DY[j]*(P[i][j]-P[i][j-1]);

            d = VIS[i][j]*((V[i][j+1]-2*V[i][j]+V[i][j-1])/(DYV[j+1]*DYV[j])+(V[i+1][j]-2*V[i][j]+V[i-1][j])/(DXV[i+1]*DXV[i]));

            t = VIS_T[i][j]*((V[i][j+1]-2*V[i][j]+V[i][j-1])/(DYV[j+1]*DYV[j])+(V[i+1][j]-2*V[i][j]+V[i-1][j])/(DXV[i+1]*DXV[i]));

            g = gy;

            rho = 0.5*(RHO[i][j-1]+RHO[i][j]);

            V_NEW[i][j] = V[i][j] + 1/rho*dts*(c1 + c2 + p + d + g + t) + V[i][j]*(2*random()-1)*I_y[i][j];
            V[i][j] = V_NEW[i][j];
        }
    }
}

void CALC_T_t( double **T_NEW, double **T, double **U, double **V, int NPX, int NPY, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV, double dts, double **RHO, double **VIS, double **CP, double **K )
{
    int i, j;
    double c1, c2, k, phi_1, phi_2;

    for(i=1;i<=NPX-2;i++)
    {
        for(j=1;j<=NPY-2;j++)
        {
            c1 = -U[i][j]*((T[i][j]-T[i-1][j])/DX[i]);
            c2 = -V[i][j]*((T[i][j]-T[i][j-1])/DY[j]);
            k = K[i][j]/(RHO[i][j]*CP[i][j])*((T[i+1][j]-2*T[i][j]+T[i-1][j])/(DX[i+1]*DX[i])+(T[i][j+1]-2*T[i][j]+T[i][j-1])/(DY[j+1]*DY[j]));
            phi_1 = VIS[i][j]/(RHO[i][j]*CP[i][j])*(2*(((U[i][j]-U[i-1][j])/DXU[i])*((U[i][j]-U[i-1][j])/DXU[i])+((V[i][j]-V[i][j-1])/DXU[i])*((V[i][j]-V[i][j-1])/DYV[j])));
            phi_2 = VIS[i][j]/(RHO[i][j]*CP[i][j])*(((V[i][j]-V[i-1][j])/DXV[i]+(U[i][j]-U[i][j-1])/DYU[j])*((V[i][j]-V[i-1][j])/DXV[i]+(U[i][j]-U[i][j-1])/DYU[j]));

            T_NEW[i][j] = T[i][j] + dts*(c1 + c2 + k + phi_1 + phi_2);
            T[i][j] = T_NEW[i][j];
        }
    }
}

