#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void GRID_BFS( double XD, double YD, int NPX, int NPY, double *X, double *Y, double *XU, double *XV, double *YU, double *YV, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV)
{
    int i, j, k;

    double rx = 1.02;
    double ry = 0.98;

    double a_ix = (XD*(1-rx))/(1-pow(rx,NPX-1));
    double a_iy = (YD/2*(1-ry))/(1-pow(ry,(NPY-1)/2));

    // X-POS CALCULATIONS

    X[0]=0;
    DX[0]=0;

    for(i=1;i<=NPX-1;i++)
    {
        DX[i]=a_ix*pow(rx,(i-1));
        X[i]=X[i-1]+DX[i];
    }

    for(i=1;i<=NPX-1;i++)
    {
        XU[i]=X[i-1]+(X[i]-X[i-1])/2;
        XV[i]=X[i];
    }

    XU[0]=-XU[1];
    XV[0]=X[0];
    XU[NPX]=X[NPX-1]+(X[NPX-1]-XU[NPX-1]);
    XV[NPX]=X[NPX-1]+(X[NPX-1]-XV[NPX-1]);

    // DXU AND DXV CALCULATIONS

    DXU[0]=0;
    DXV[0]=0;

    for(i=1;i<=NPX;i++)
    {
        DXU[i]=XU[i]-XU[i-1];
        DXV[i]=DX[i];
    }

    // Y-POS CALCULATIONS

    Y[0]=0;
    DY[0]=0;

    for(j=1;j<=(NPY-1)/2;j++)
    {
        DY[j]=a_iy*pow(ry,(j-1));
        Y[j]=Y[j-1]+DY[j];
    }

    k = 0;

    for(j=(NPY-1)/2+1;j<=NPY-1;j++)
    {
        DY[j]=DY[j-1-2*k];
        Y[j]=Y[j-1]+DY[j];
        k = k + 1;
    }

    for(j=1;j<=NPY-1;j++)
    {
        YU[j]=Y[j];
        YV[j]=Y[j-1]+(Y[j]-Y[j-1])/2;
    }

    YU[0]=Y[0];
    YV[0]=-YV[1];
    YU[NPY]=Y[NPY-1]+(Y[NPY-1]-YU[NPY-1]);
    YV[NPY]=Y[NPY-1]+(Y[NPY-1]-YV[NPY-1]);

    // DYU AND DYV CALCULATIONS

    DYU[0]=0;
    DYV[0]=0;

    for(j=1;j<=NPY;j++)
    {
        DYU[j]=DY[j];
        DYV[j]=YV[j]-YV[j-1];
    }
}

void GRID_2( double XD, double YD, int NPX, int NPY, double *X, double *Y, double *XU, double *XV, double *YU, double *YV, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV)
{
    int i, j;

    double rx = 1.0001;
    double ry = 1.0100;

    double a_ix = (XD*(1-rx))/(1-pow(rx,NPX-1));
    double a_iy = (YD*(1-ry))/(1-pow(ry,NPY-1));

    // X-POS CALCULATIONS

    X[0]=0;
    DX[0]=0;

    for(i=1;i<=NPX-1;i++)
    {
        DX[i]=a_ix*pow(rx,(i-1));
        X[i]=X[i-1]+DX[i];
    }

    for(i=1;i<=NPX-1;i++)
    {
        XU[i]=X[i-1]+(X[i]-X[i-1])/2;
        XV[i]=X[i];
    }

    XU[0]=-XU[1];
    XV[0]=X[0];
    XU[NPX]=X[NPX-1]+(X[NPX-1]-XU[NPX-1]);
    XV[NPX]=X[NPX-1]+(X[NPX-1]-XV[NPX-1]);

    // DXU AND DXV CALCULATIONS

    DXU[0]=0;
    DXV[0]=0;

    for(i=1;i<=NPX;i++)
    {
        DXU[i]=XU[i]-XU[i-1];
        DXV[i]=DX[i];
    }

    // Y-POS CALCULATIONS

    Y[0]=0;
    DY[0]=0;

    for(j=1;j<=NPY-1;j++)
    {
        DY[j]=a_iy*pow(ry,(j-1));
        Y[j]=Y[j-1]+DY[j];
    }

    for(j=1;j<=NPY-1;j++)
    {
        YU[j]=Y[j];
        YV[j]=Y[j-1]+(Y[j]-Y[j-1])/2;
    }

    YU[0]=Y[0];
    YV[0]=-YV[1];
    YU[NPY]=Y[NPY-1]+(Y[NPY-1]-YU[NPY-1]);
    YV[NPY]=Y[NPY-1]+(Y[NPY-1]-YV[NPY-1]);

    // DYU AND DYV CALCULATIONS

    DYU[0]=0;
    DYV[0]=0;

    for(j=1;j<=NPY;j++)
    {
        DYU[j]=DY[j];
        DYV[j]=YV[j]-YV[j-1];
    }
}

void GRID_3( double XD, double YD, int NPX, int NPY, double *X, double *Y, double *XU, double *XV, double *YU, double *YV, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV)
{
    int i, j, k;

    double rx = 0.975;
    double ry = 1.025;

    double a_ix = (XD/2*(1-rx))/(1-pow(rx,(NPX-1)/2));
    double a_iy = (YD*(1-ry))/(1-pow(ry,NPY-1));

    // X-POS CALCULATIONS

    X[0]=0;
    DX[0]=0;

    for(i=1;i<=(NPX-1)/2;i++)
    {
        DX[i]=a_ix*pow(rx,(i-1));
        X[i]=X[i-1]+DX[i];
    }

    k = 0;

    for(i=(NPX-1)/2+1;i<=NPX-1;i++)
    {
        DX[i]=DX[i-1-2*k];
        X[i]=X[i-1]+DX[i];
        k = k + 1;
    }

    for(i=1;i<=NPX-1;i++)
    {
        XU[i]=X[i-1]+(X[i]-X[i-1])/2;
        XV[i]=X[i];
    }

    XU[0]=-XU[1];
    XV[0]=X[0];
    XU[NPX]=X[NPX-1]+(X[NPX-1]-XU[NPX-1]);
    XV[NPX]=X[NPX-1]+(X[NPX-1]-XV[NPX-1]);

    // DXU AND DXV CALCULATIONS

    DXU[0]=0;
    DXV[0]=0;

    for(i=1;i<=NPX;i++)
    {
        DXU[i]=XU[i]-XU[i-1];
        DXV[i]=DX[i];
    }

    // Y-POS CALCULATIONS

    Y[0]=0;
    DY[0]=0;

    for(j=1;j<=NPY-1;j++)
    {
        DY[j]=a_iy*pow(ry,(j-1));
        Y[j]=Y[j-1]+DY[j];
    }

    for(j=1;j<=NPY-1;j++)
    {
        YU[j]=Y[j];
        YV[j]=Y[j-1]+(Y[j]-Y[j-1])/2;
    }

    YU[0]=Y[0];
    YV[0]=-YV[1];
    YU[NPY]=Y[NPY-1]+(Y[NPY-1]-YU[NPY-1]);
    YV[NPY]=Y[NPY-1]+(Y[NPY-1]-YV[NPY-1]);

    // DYU AND DYV CALCULATIONS

    DYU[0]=0;
    DYV[0]=0;

    for(j=1;j<=NPY;j++)
    {
        DYU[j]=DY[j];
        DYV[j]=YV[j]-YV[j-1];
    }
}

void GRID_CAVITY( double XD, double YD, int NPX, int NPY, double *X, double *Y, double *XU, double *XV, double *YU, double *YV, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV)
{
    int i, j, k;

    double rx = 1.01;
    double ry = 1.01;

    double a_ix = (XD/2*(1-rx))/(1-pow(rx,(NPX-1)/2));
    double a_iy = (YD/2*(1-ry))/(1-pow(ry,(NPY-1)/2));

    // X-POS CALCULATIONS

    X[0]=0;
    DX[0]=0;

    for(i=1;i<=(NPX-1)/2;i++)
    {
        DX[i]=a_ix*pow(rx,(i-1));
        X[i]=X[i-1]+DX[i];
    }

    k = 0;

    for(i=(NPX-1)/2+1;i<=NPX-1;i++)
    {
        DX[i]=DX[i-1-2*k];
        X[i]=X[i-1]+DX[i];
        k = k + 1;
    }

    for(i=1;i<=NPX-1;i++)
    {
        XU[i]=X[i-1]+(X[i]-X[i-1])/2;
        XV[i]=X[i];
    }

    XU[0]=-XU[1];
    XV[0]=X[0];
    XU[NPX]=X[NPX-1]+(X[NPX-1]-XU[NPX-1]);
    XV[NPX]=X[NPX-1]+(X[NPX-1]-XV[NPX-1]);

    // DXU AND DXV CALCULATIONS

    DXU[0]=0;
    DXV[0]=0;

    for(i=1;i<=NPX;i++)
    {
        DXU[i]=XU[i]-XU[i-1];
        DXV[i]=DX[i];
    }

    // Y-POS CALCULATIONS

    Y[0]=0;
    DY[0]=0;

    for(j=1;j<=(NPY-1)/2;j++)
    {
        DY[j]=a_iy*pow(ry,(j-1));
        Y[j]=Y[j-1]+DY[j];
    }

    k = 0;

    for(j=(NPY-1)/2+1;j<=NPY-1;j++)
    {
        DY[j]=DY[j-1-2*k];
        Y[j]=Y[j-1]+DY[j];
        k = k + 1;
    }

    for(j=1;j<=NPY-1;j++)
    {
        YU[j]=Y[j];
        YV[j]=Y[j-1]+(Y[j]-Y[j-1])/2;
    }

    YU[0]=Y[0];
    YV[0]=-YV[1];
    YU[NPY]=Y[NPY-1]+(Y[NPY-1]-YU[NPY-1]);
    YV[NPY]=Y[NPY-1]+(Y[NPY-1]-YV[NPY-1]);

    // DYU AND DYV CALCULATIONS

    DYU[0]=0;
    DYV[0]=0;

    for(j=1;j<=NPY;j++)
    {
        DYU[j]=DY[j];
        DYV[j]=YV[j]-YV[j-1];
    }
}

void GRID( double XD, double YD, int NPX, int NPY, double *X, double *Y, double *XU, double *XV, double *YU, double *YV, double *DX, double *DY, double *DXU, double *DXV, double *DYU, double *DYV)
{
    int i, j, k;

    double rx = 1.00001;
    double ry = 0.98;

    double a_ix = (XD*(1-rx))/(1-pow(rx,NPX-1));
    double a_iy = (YD/2*(1-ry))/(1-pow(ry,(NPY-1)/2));

    // X-POS CALCULATIONS

    X[0]=0;
    DX[0]=0;

    for(i=1;i<=NPX-1;i++)
    {
        DX[i]=a_ix*pow(rx,(i-1));
        X[i]=X[i-1]+DX[i];
    }

    for(i=1;i<=NPX-1;i++)
    {
        XU[i]=X[i-1]+(X[i]-X[i-1])/2;
        XV[i]=X[i];
    }

    XU[0]=-XU[1];
    XV[0]=X[0];
    XU[NPX]=X[NPX-1]+(X[NPX-1]-XU[NPX-1]);
    XV[NPX]=X[NPX-1]+(X[NPX-1]-XV[NPX-1]);

    // DXU AND DXV CALCULATIONS

    DXU[0]=0;
    DXV[0]=0;

    for(i=1;i<=NPX;i++)
    {
        DXU[i]=XU[i]-XU[i-1];
        DXV[i]=DX[i];
    }

    // Y-POS CALCULATIONS

    Y[0]=0;
    DY[0]=0;

    for(j=1;j<=(NPY-1)/2;j++)
    {
        DY[j]=a_iy*pow(ry,(j-1));
        Y[j]=Y[j-1]+DY[j];
    }

    k = 0;

    for(j=(NPY-1)/2+1;j<=NPY-1;j++)
    {
        DY[j]=DY[j-1-2*k];
        Y[j]=Y[j-1]+DY[j];
        k = k + 1;
    }

    for(j=1;j<=NPY-1;j++)
    {
        YU[j]=Y[j];
        YV[j]=Y[j-1]+(Y[j]-Y[j-1])/2;
    }

    YU[0]=Y[0];
    YV[0]=-YV[1];
    YU[NPY]=Y[NPY-1]+(Y[NPY-1]-YU[NPY-1]);
    YV[NPY]=Y[NPY-1]+(Y[NPY-1]-YV[NPY-1]);

    // DYU AND DYV CALCULATIONS

    DYU[0]=0;
    DYV[0]=0;

    for(j=1;j<=NPY;j++)
    {
        DYU[j]=DY[j];
        DYV[j]=YV[j]-YV[j-1];
    }
}
