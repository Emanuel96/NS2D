#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void BC_BFS(double **U, double **V, double **P, int NPX, int NPY, double U_in, double *X, double *Y, double time, double **T)
{
    int i, j;
    double Y_STEP = 0;
    //double a = PID_c;

    //CALC_U

    for(j=0;j<=NPY-1;j++)
    {
        U[1][j]=0;
        //U[1][j]=U[2][j];
        U[0][j]=U[1][j];

        if ( Y[j] >= Y_STEP && Y[j] <= Y[NPY-1]-Y_STEP )
        {
            U[1][j]=U_in;
            //U[1][j]=U[2][j];
            //U[1][j]=U_in+U_in*sin(2*3.1415*(1/0.5)*3.5*Y[j]);
            //U[1][j]=(U_in/log(a*0.5+1))*log(a*(Y[j]-0.00)+1);
            U[0][j]=U[1][j];
        }
    }

    for(j=0;j<=NPY-1;j++)
    {
        U[NPX-1][j]=U[NPX-2][j];
        U[NPX][j]=U[NPX-1][j];
    }

    for(i=0;i<=NPX-1;i++)
    {
        U[i][0]=0;
        U[i][NPY-1]=0;//U[i][NPY-2];//0;
    }

    //CALC_V

    for(j=0;j<=NPY-1;j++)
    {
        V[0][j]=0;//V[1][j];//0;
        V[NPX-1][j]=V[NPX-2][j];
    }

    //CALC_P

    for(i=0;i<=NPX-1;i++)
    {
        P[i][0]=P[i][1];
        P[i][NPY-1]=P[i][NPY-2];
    }

    for(j=0;j<=NPY-1;j++)
    {
        P[0][j]=P[1][j];
        P[NPX-1][j]=P[NPX-2][j];
    }

    //CALC_T

//    for(i=0;i<=NPX-1;i++)
//    {
//        for(j=0;j<=NPY-1;j++)
//        {
//            if ( ( ( X[i] >= 0.25 && X[i] <= 0.50 ) && ( Y[j] >= 0.4 && Y[j] <= 0.6 ) ) )
//            {
//                T[i][j] = 273.15 - 100/(1+exp(-0.5*(time-10)));
//            }
//        }
//
//        //T[i][0]=50;
//        //T[i][NPY-1]=300;
//    }
//
//    for(j=0;j<=NPY-1;j++)
//    {
//        //T[0][j] = 300;
//        T[NPX-1][j]=T[NPX-2][j];
//    }

}

void BC(double **U, double **V, double **P, int NPX, int NPY, double U_in, double *X, double *Y, double time, double **T, double PID_c)
{
    int i, j;
    double Y_STEP = 0;
    //double a = PID_c;

    //CALC_U

    for(j=0;j<=NPY-1;j++)
    {
        U[1][j]=0;
        //U[1][j]=U[2][j];
        U[0][j]=U[1][j];

        if ( Y[j] >= Y_STEP && Y[j] <= Y[NPY-1]-Y_STEP )
        {
            U[1][j]=U_in;
            //U[1][j]=U[2][j];
            //U[1][j]=U_in+U_in*sin(2*3.1415*(1/0.5)*3.5*Y[j]);
            //U[1][j]=(U_in/log(a*0.5+1))*log(a*(Y[j]-0.00)+1);
            U[0][j]=U[1][j];
        }
    }

    for(j=0;j<=NPY-1;j++)
    {
        U[NPX-1][j]=U[NPX-2][j];
        U[NPX][j]=U[NPX-1][j];
    }

    for(i=0;i<=NPX-1;i++)
    {
        U[i][0]=0;
        U[i][NPY-1]=0;//U[i][NPY-2];//0;
    }

    //CALC_V

    for(j=0;j<=NPY-1;j++)
    {
        V[0][j]=0;//V[1][j];//0;
        V[NPX-1][j]=V[NPX-2][j];
    }

    //CALC_P

    for(i=0;i<=NPX-1;i++)
    {
        P[i][0]=P[i][1];
        P[i][NPY-1]=P[i][NPY-2];
    }

    for(j=0;j<=NPY-1;j++)
    {
        P[0][j]=P[1][j];
        P[NPX-1][j]=P[NPX-2][j];
    }

    //CALC_T

//    for(i=0;i<=NPX-1;i++)
//    {
//        T[i][0]=300;
//        T[i][NPY-1]=300;
//    }

//    for(j=0;j<=NPY-1;j++)
//    {
//        T[0][j] = 300;
//        T[NPX-1][j]=T[NPX-2][j];
//    }

//    for(i=0;i<=NPX-1;i++)
//    {
//        for(j=0;j<=NPY-1;j++)
//        {
//            if( Y[j] <= 0.2*exp(-20*pow(X[i]-0.75,2))*1/(1+exp(-3*(time-3))) || Y[j] >= 1-(0.2*exp(-20*pow(X[i]-0.75,2)))*1/(1+exp(-3*(time-3))) || ( X[i] >= 0.75 && Y[j] >= 1-0.2*(1/(1+exp(-3*(time-3)))) ) || ( X[i] >= 0.75 && Y[j] <= 0.2*(1/(1+exp(-3*(time-3)))) ) )
//            //if( Y[j] <= 0.3*exp(-20*pow(X[i]-0.75,2))*1/(1+exp(-3*(time-3))) || Y[j] >= 1-(0.3*exp(-20*pow(X[i]-0.75,2)))*1/(1+exp(-3*(time-3))) || ( X[i] >= 0.75 && Y[j] >= 1-0.3*(1/(1+exp(-3*(time-3)))) ) || ( X[i] >= 0.75 && Y[j] <= 0.3*(1/(1+exp(-3*(time-3)))) ) )
//            //if( Y[j] >= 1-0.1*exp(-20*pow(X[i]-0.75,2)) )//*sin(2*3.1415*1/20*time*0.25) )// || Y[j] >= 0.5-0.15/(1+exp(100000*(X[i]-0.25)))*sin(2*3.1415*1/20*time*0.25) )
//            {
//                U[i][j]=0;
//                U[i+1][j]=0;
//                V[i][j]=0;
//                V[i][j+1]=0;
//            }
//        }
//    }

    //for(i=0;i<=NPX-1;i++)
    //{
        //T[i][0]=(500/log(a*1.5+1))*log(a*X[i]+1);
        //T[i][0]=500*sin(2*3.1415*(1/1.5)*5*X[i]);
    //}

    // TESTING

//    double a_x, a_y, f;
//    a_x = 0.00;
//    a_y = 0.05;
//    f = 0.25;
//
//    for(i=0;i<=NPX-1;i++)
//    {
//        for(j=0;j<=NPY-1;j++)
//        {
//            //if ( ( ( X[i] >= 0.5 && X[i] <= 1.0 ) && ( Y[j] >= 0.4 && Y[j] <= 0.6 ) ) )// || ( (X[i] >= 0.0 && X[i] <= 0.2) && (Y[j] >= 0.0 && Y[j] <= 0.05) ) )
//            if ( pow(X[i]-1+a_x*sin(2*3.1415*f*time),2)/(0.5*0.5) + pow(Y[j]-0.5+a_y*cos(2*3.1415*f*time),2)/(0.05*0.05) <= 1 )
//            {
//                U[i][j]=2*3.1415*f*a_x*cos(2*3.1415*f*time);
//                U[i+1][j]=2*3.1415*f*a_x*cos(2*3.1415*f*time);
//                V[i][j]=-2*3.1415*f*a_y*sin(2*3.1415*f*time);
//                V[i][j+1]=-2*3.1415*f*a_y*sin(2*3.1415*f*time);
//            }
//        }
//    }

//    for(i=0;i<=NPX-1;i++)
//    {
//        for(j=0;j<=NPY-1;j++)
//        {
//            if ( ( ( X[i] >= 0.5 && X[i] <= 1.0 ) && ( Y[j] >= 0.1 && Y[j] <= 0.4 ) ) )// || ( (X[i] >= 0.0 && X[i] <= 0.2) && (Y[j] >= 0.0 && Y[j] <= 0.05) ) )
//            //if ( pow(X[i]-0.75,2)/(0.25*0.25) + pow(Y[j]+0.05*sin(2*3.1415*0.1*time)-0.5,2)/(0.1*0.1) <= 1 )
//            {
//                T[i][j] = 500;
//                //U[i][j]=0;
//                //U[i+1][j]=0;
//                //V[i][j]=2*3.1415*0.1*0.05*cos(2*3.1415*0.1*time);
//                //V[i][j+1]=2*3.1415*0.1*0.05*cos(2*3.1415*0.1*time);
//            }
//        }
//    }

//    for(i=0;i<=NPX-1;i++)
//    {
//        for(j=0;j<=NPY-1;j++)
//        {
//            if ( Y[j] <= (0.10*1*PID_c)/(1+exp(-10*(X[i]-0.75))) )// || Y[j] >= 0.5-(0.10*PID_c)/(1+exp(-10*(X[i]-0.75))) )// || Y[j] <= 0.15/(1+exp(-100000*(X[i]-0.5))) )//*sin(2*3.1415*1/20*time*0.25) )// || Y[j] >= 0.5-0.15/(1+exp(100000*(X[i]-0.25)))*sin(2*3.1415*1/20*time*0.25) )
//            //if( Y[j] <= PID_c*X[i] || Y[j] >= 0.5-PID_c*X[i] )
//            {
//                U[i][j]=0;
//                U[i+1][j]=0;
//                V[i][j]=0;
//                V[i][j+1]=0;
//                //T[i][j]=300;
//            }
//        }
//    }
}

void BC_1(double **U, double **V, double **P, int NPX, int NPY, double U_in, double *X, double *Y, double time, double **T, double PID_c)
{
    int i, j;

    //CALC_U

    for(j=0;j<=NPY-1;j++)
    {
        //U[1][j]=((0.5*U_in/log(500*0.5+1))*log(500*(Y[j]-0.00)+1))*exp(-10*pow(time-5,2));
        U[1][j]=U[2][j];
        U[0][j]=U[1][j];
        U[NPX-1][j]=U[NPX-2][j];
        U[NPX][j]=U[NPX-1][j];
    }

    for(i=0;i<=NPX-1;i++)
    {
        U[i][0]=0;
        U[i][NPY-1]=0;
    }

    //CALC_V

    double jet_pos = 0.50;
    double delta_jet = 0.05;

    for(i=0;i<=NPX-1;i++)
    {
        if ( X[i] >= jet_pos-delta_jet && X[i] <= jet_pos+delta_jet )
        {
            V[i][NPY]=-U_in;
            V[i][NPY-1]=-U_in;
            //T[i][NPY] = 400;
        }
    }

    for(j=0;j<=NPY-1;j++)
    {
        V[0][j]=V[i][j];
        V[NPX-1][j]=V[NPX-2][j];
    }

    //CALC_P

    for(i=0;i<=NPX-1;i++)
    {
        P[i][0]=P[i][1];
        P[i][NPY-1]=P[i][NPY-2];
    }

    for(j=0;j<=NPY-1;j++)
    {
        P[0][j]=P[1][j];
        P[NPX-1][j]=P[NPX-2][j];
    }
}

void BC_SIN(double **U, double **V, double **P, int NPX, int NPY, double U_in, double *X, double *Y, double time, double **T)
{
    int i, j;
    double Y_STEP = 0.0;

    //CALC_U

    for(j=1;j<=NPY-2;j++)
    {
        U[1][j]=0;
        //U[1][j]=U[2][j];
        U[0][j]=U[1][j];

        if ( Y[j] >= Y_STEP && Y[j] <= Y[NPY-1]-Y_STEP )
        {
            U[1][j]=U_in;
            //U[1][j]=U[2][j];
            U[0][j]=U[1][j];
        }
    }

    for(j=0;j<=NPY-1;j++)
    {
        U[NPX-1][j]=U[NPX-2][j];
        U[NPX][j]=U[NPX-1][j];
    }

    for(i=0;i<=NPX-1;i++)
    {
        U[i][0]=0;
        U[i][NPY-1]=U[i][NPY-2];
    }

    //CALC_V

    for(j=0;j<=NPY-1;j++)
    {
        V[0][j]=0;
        V[NPX-1][j]=V[NPX-2][j];
    }

    // TESTING

    for(i=0;i<=NPX-1;i++)
    {
        for(j=0;j<=NPY-1;j++)
        {
            if ( Y[j] <= 0.1*sin(2*3.1415*1.5/4*X[i]) || Y[j] >= 0.5-0.1*sin(2*3.1415*1.5/4*X[i]) )
            {
                U[i][j]=0;
                U[i+1][j]=0;
                V[i][j]=0;
                V[i][j+1]=0;
            }
        }
    }

    //CALC_P

    for(i=0;i<=NPX-1;i++)
    {
        P[i][0]=P[i][1];
        P[i][NPY-1]=P[i][NPY-2];
    }

    for(j=0;j<=NPY-1;j++)
    {
        P[0][j]=P[1][j];
        P[NPX-1][j]=P[NPX-2][j];
    }
}

void BC_DONTKNOW(double **U, double **V, double **P, int NPX, int NPY, double U_in, double *X, double *Y, double time, double **T)
{
    int i, j;
    double Y_STEP=0.0;

    //CALC_U

    for(j=1;j<=NPY-2;j++)
    {
        U[1][j]=0;
        //U[1][j]=U[2][j];
        U[0][j]=U[1][j];

        if ( Y[j] >= Y_STEP && Y[j] <= Y[NPY-1]-Y_STEP )
        {
            U[1][j]=U_in;
            //U[1][j]=U[2][j];
            U[0][j]=U[1][j];
        }
    }

    for(j=0;j<=NPY-1;j++)
    {
        U[NPX-1][j]=U[NPX-2][j];
        U[NPX][j]=U[NPX-1][j];
    }

    for(i=0;i<=NPX-1;i++)
    {
        U[i][0]=0;
        U[i][NPY-1]=0;
    }

    //CALC_V

    for(j=0;j<=NPY-1;j++)
    {
        V[0][j]=0;
        V[NPX-1][j]=V[NPX-2][j];
    }

    //CALC_P

    for(i=0;i<=NPX-1;i++)
    {
        P[i][0]=P[i][1];
        P[i][NPY-1]=P[i][NPY-2];
    }

    for(j=0;j<=NPY-1;j++)
    {
        P[0][j]=P[1][j];
        P[NPX-1][j]=P[NPX-2][j];
    }

    for(i=0;i<=NPX-1;i++)
    {
        for(j=0;j<=NPY-1;j++)
        {
            //if ( ( (X[i] >= 0.3 && X[i] <= 1.0) && (Y[j] >= 0.0 && Y[j] <= 0.05) ) || ( (X[i] >= 0.0 && X[i] <= 0.2) && (Y[j] >= 0.0 && Y[j] <= 0.05) ) )
            if ( (X[i]-0.2)*(X[i]-0.2)/(0.1*0.1) + (Y[j]+0.025*sin(2*3.1415*0.5*time)*(80*X[i]*X[i]*X[i]-0.125)-0.125)*(Y[j]+0.025*sin(2*3.1415*0.5*time)*(80*X[i]*X[i]*X[i]-0.125)-0.125)/(0.01*0.01) <= 1 )
            {
                U[i][j]=0;
                U[i+1][j]=0;
                V[i][j]=-2*3.1415*0.01*0.5*cos(2*3.1415*0.5*time);
                V[i][j+1]=-2*3.1415*0.01*0.5*cos(2*3.1415*0.5*time);
            }
        }
    }
}

void BC_CAVITY(double **U, double **V, double **P, int NPX, int NPY, double U_in, double *X, double *Y, double time, double **T)
{
    int i, j;

    //CALC_U

    for(i=0;i<=NPX;i++)
    {
        U[i][0]=U_in;//0;
        U[i][NPY-1]=0;
    }

    //CALC_V

    for(j=0;j<=NPY;j++)
    {
        V[0][j]=0;
        V[NPX-1][j]=0;
    }

    //CALC_P

    for(i=0;i<=NPX-1;i++)
    {
        P[i][0]=P[i][1];
        P[i][NPY-1]=P[i][NPY-2];
    }

    for(j=0;j<=NPY-1;j++)
    {
        P[0][j]=P[1][j];
        P[NPX-1][j]=P[NPX-2][j];
    }

    //CALC_T

    for(i=0;i<=NPX-1;i++)
    {
        T[i][0]=500;
    }
}


void BC_SBY(double **U, double **V, double **P, int NPX, int NPY, double U_in, double *X, double *Y, double time, double **T)
{
    int i, j;
    double Y_STEP=0.0;

    //CALC_U

    for(j=1;j<=NPY-2;j++)
    {
        U[1][j]=0;
        //U[1][j]=U[2][j];
        U[0][j]=U[1][j];

        if ( Y[j] >= Y_STEP && Y[j] <= Y[NPY-1]-Y_STEP )
        {
            U[1][j]=U_in;
            //U[1][j]=U[2][j];
            U[0][j]=U[1][j];
        }
    }

    for(j=0;j<=NPY-1;j++)
    {
        U[NPX-1][j]=U[NPX-2][j];
        U[NPX][j]=U[NPX-1][j];
    }

    for(i=0;i<=NPX-1;i++)
    {
        U[i][0]=0;
        U[i][NPY-1]=0;
    }

    //CALC_V

    for(j=0;j<=NPY-1;j++)
    {
        V[0][j]=0;
        V[NPX-1][j]=V[NPX-2][j];
    }

//    for(i=0;i<=NPX-1;i++)
//    {
//        if ( X[i] >= 0.225 && X[i] <= 0.275)
//        {
//            V[i][NPY-1]=-U_in;
//            V[i][NPY]=V[i][NPY-1];
//        }
//    }

    //CALC_P

    for(i=0;i<=NPX-1;i++)
    {
        P[i][0]=P[i][1];
        P[i][NPY-1]=P[i][NPY-2];
    }

    for(j=0;j<=NPY-1;j++)
    {
        P[0][j]=P[1][j];
        P[NPX-1][j]=P[NPX-2][j];
    }

    //CALC_T

//    for(i=0;i<=NPX-1;i++)
//    {
//        if (X[i] >= 0.05 && X[i] <= 0.10)
//        {
//            //U[i][50]=0;
//            //U[i+1][50]=0;
//            //V[i][50]=0;
//            //V[i][50+1]=0;
//            //T[i][50]=500;
//        }
//    }

    for(j=0;j<=NPY-1;j++)
    {
        //T[0][j]=150;
        //T[0][j]=150*exp(-0.01*time);
        T[NPX-1][j]=T[NPX-2][j];
    }

//    for(i=0;i<=NPX-1;i++)
//    {
//        for(j=0;j<=NPY-1;j++)
//        {
//            if ( (X[i] >= 0.15 && X[i] <= 0.25 ) && (Y[j] >= 0.1 && Y[j] <= 0.15) )
//            if ( ( (X[i] >= 0.3 && X[i] <= 1.0) && (Y[j] >= 0.0 && Y[j] <= 0.05) ) || ( (X[i] >= 0.0 && X[i] <= 0.2) && (Y[j] >= 0.0 && Y[j] <= 0.05) ) )
//            if ( (X[i]-0.2)*(X[i]-0.2) + (Y[j]-0.125)*(Y[j]-0.125) <= 0.05*0.05 )
//            {
//                U[i][j]=0;
//                U[i+1][j]=0;
//                V[i][j]=0;
//                V[i][j+1]=0;
//                T[i][j]=1000;
//            }
//        }
//    }

    for(i=0;i<=NPX-1;i++)
    {
        for(j=0;j<=NPY-1;j++)
        {
            //if ( ( (X[i] >= 0.3 && X[i] <= 1.0) && (Y[j] >= 0.0 && Y[j] <= 0.05) ) || ( (X[i] >= 0.0 && X[i] <= 0.2) && (Y[j] >= 0.0 && Y[j] <= 0.05) ) )
//            if ( (X[i]-0.2)*(X[i]-0.2)/(0.05*0.05) + (Y[j]+0.01*sin(2*3.1415*0.5*time)-0.125)*(Y[j]+0.01*sin(2*3.1415*0.5*time)-0.125)/(0.05*0.05) <= 1 )
//            {
//                U[i][j]=0;
//                U[i+1][j]=0;
//                V[i][j]=-2*3.1415*0.01*0.5*cos(2*3.1415*0.5*time);
//                V[i][j+1]=-2*3.1415*0.01*0.5*cos(2*3.1415*0.5*time);
//                //T[i][j]=500;
//            }
//            if ( (X[i]-0.2)*(X[i]-0.2)/(0.025*0.025) + (Y[j]-0.125)*(Y[j]-0.125)/(0.05*0.05) <= 1 )
//            {
//                U[i][j]=0;
//                U[i+1][j]=0;
//                V[i][j]=0;
//                V[i][j+1]=0;
//                //T[i][j]=500;
//            }
        }
    }

}

void BC_CYLINDER(double **U, double **V, double **P, int NPX, int NPY, double U_in, double *X, double *Y, double time, double **T)
{
    int i, j;
    double Y_STEP = 0;
    //double a = PID_c;

    //CALC_U

    for(j=0;j<=NPY-1;j++)
    {
        U[1][j]=0;
        //U[1][j]=U[2][j];
        U[0][j]=U[1][j];

        if ( Y[j] >= Y_STEP && Y[j] <= Y[NPY-1]-Y_STEP )
        {
            U[1][j]=U_in;
            //U[1][j]=U[2][j];
            //U[1][j]=U_in+U_in*sin(2*3.1415*(1/0.5)*3.5*Y[j]);
            //U[1][j]=(U_in/log(a*0.5+1))*log(a*(Y[j]-0.00)+1);
            U[0][j]=U[1][j];
        }
    }

    for(j=0;j<=NPY-1;j++)
    {
        U[NPX-1][j]=0;//U[NPX-2][j];
        U[NPX][j]=0;//U[NPX-1][j];
    }

    for(i=0;i<=NPX-1;i++)
    {
        U[i][0]=0;
        U[i][NPY-1]=0;//U[i][NPY-2];//0;
    }

    //CALC_V

    for(j=0;j<=NPY-1;j++)
    {
        V[0][j]=0;//V[1][j];//0;
        V[NPX-1][j]=0;//V[NPX-2][j];
    }

    //CALC_P

    for(i=0;i<=NPX-1;i++)
    {
        P[i][0]=P[i][1];
        P[i][NPY-1]=P[i][NPY-2];
    }

    for(j=0;j<=NPY-1;j++)
    {
        P[0][j]=P[1][j];
        P[NPX-1][j]=P[NPX-2][j];
    }

    //CYLINDER

    for(i=0;i<=NPX-1;i++)
    {
        for(j=0;j<=NPY-1;j++)
        {
            if ( (X[i]+0.05*cos(2*3.1415*0.5*time)-0.25)*(X[i]+0.05*cos(2*3.1415*0.5*time)-0.25)/(0.05*0.05) + (Y[j]+0.05*sin(2*3.1415*0.5*time)-0.25)*(Y[j]+0.0*sin(2*3.1415*0.5*time)-0.25)/(0.05*0.05) <= 1 )
            {
                U[i][j]=2*3.1415*0.05*0.5*sin(2*3.1415*0.5*time);
                U[i+1][j]=U[i][j];
                V[i][j]=-2*3.1415*0.05*0.5*cos(2*3.1415*0.5*time);
                V[i][j+1]=V[i][j];
                //T[i][j]=1000;
            }
        }
    }

}
