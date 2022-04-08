#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "functions.h"
#include "INIT.h"
#include "GRID.h"
#include "BCs.h"
#include "SOLVER.h"
#include "TDMA.h"
#include "CALC_P.h"
#include "PROPS.h"
#include "POSPRO.h"
#include "MDO.h"
#include <string.h>

#define NL 1000
#define NC 1000

int main()
{
    //Variables

    FILE* file;
    FILE* file_R;
    FILE* file_MDO;
    FILE* file_in;

    double **U = (double **) malloc(NL * sizeof(double*));
    double **U_NEW = (double **) malloc(NL * sizeof(double*));
    double **V = (double **) malloc(NL * sizeof(double*));
    double **V_NEW = (double **) malloc(NL * sizeof(double*));
    double **P = (double **) malloc(NL * sizeof(double*));
    double **PC = (double **) malloc(NL * sizeof(double*));
    double **R = (double **) malloc(NL * sizeof(double*));

    double **RHO = (double **) malloc(NL * sizeof(double*));
    double **VIS = (double **) malloc(NL * sizeof(double*));
    double **VIS_T = (double **) malloc(NL * sizeof(double*));

    double **T = (double **) malloc(NL * sizeof(double*));
    double **T_NEW = (double **) malloc(NL * sizeof(double*));
    double **CP = (double **) malloc(NL * sizeof(double*));
    double **K = (double **) malloc(NL * sizeof(double*));

    double **Re = (double **) malloc(NL * sizeof(double*));

    double **I_x = (double **) malloc(NL * sizeof(double*));
    double **I_y = (double **) malloc(NL * sizeof(double*));

    for (int i=0; i<NL; i++)
    {
        U[i] = (double *)malloc(NC * sizeof(double));
        U_NEW[i] = (double *)malloc(NC * sizeof(double));
        V[i] = (double *)malloc(NC * sizeof(double));
        V_NEW[i] = (double *)malloc(NC * sizeof(double));
        P[i] = (double *)malloc(NC * sizeof(double));
        PC[i] = (double *)malloc(NC * sizeof(double));
        R[i] = (double *)malloc(NC * sizeof(double));

        RHO[i] = (double *)malloc(NC * sizeof(double));
        VIS[i] = (double *)malloc(NC * sizeof(double));
        VIS_T[i] = (double *)malloc(NC * sizeof(double));

        T[i] = (double *)malloc(NC * sizeof(double));
        T_NEW[i] = (double *)malloc(NC * sizeof(double));
        CP[i] = (double *)malloc(NC * sizeof(double));
        K[i] = (double *)malloc(NC * sizeof(double));

        Re[i] = (double *)malloc(NC * sizeof(double));

        I_x[i] = (double *)malloc(NC * sizeof(double));
        I_y[i] = (double *)malloc(NC * sizeof(double));
    }

    double *X = (double *) malloc(NC * sizeof(double*));
    double *Y = (double *) malloc(NL * sizeof(double*));
    double *XU = (double *) malloc(NC * sizeof(double*));
    double *XV = (double *) malloc(NC * sizeof(double*));
    double *YU = (double *) malloc(NL * sizeof(double*));
    double *YV = (double *) malloc(NL * sizeof(double*));
    double *DX = (double *) malloc(NC * sizeof(double*));
    double *DY = (double *) malloc(NL * sizeof(double*));
    double *DXU = (double *) malloc(NC * sizeof(double*));
    double *DXV = (double *) malloc(NC * sizeof(double*));
    double *DYU = (double *) malloc(NL * sizeof(double*));
    double *DYV = (double *) malloc(NL * sizeof(double*));

    //TDMA SOLVER

    double *A = (double *) malloc(NC * sizeof(double*));
    double *B = (double *) malloc(NC * sizeof(double*));
    double *D = (double *) malloc(NC * sizeof(double*));
    double *C = (double *) malloc(NC * sizeof(double*));

    double w, C1;
    double urf;
    int control;

    char fname[100];

    printf("Navier-Stokes Solver 2D V4.0.0. (2020)\n");
    printf("Emanuel Camacho. All rights reserved.\n\n");

    //strcpy(fname,"INPUT.dat");
    //strcpy(fname,"INPUT_CAVITY.dat");
    strcpy(fname,"./INPUT/INPUT_BFS.dat");

    //printf("Insert the input file name: ");
    //gets(fname);

    file_in = fopen(fname,"r");

    if( file_in == NULL )
    {
        printf("[ERROR] Input file(s) was/were not found.\n\n");
        getchar();
        return 0;
    }

    char c[100];

    void readString()
    {
        fscanf(file_in,"%[^\n]", c);
    }

    int NX, NY, NI, NPX, NPY, NTS;
    double dts, dt;

    readString();
    fscanf(file_in,"%d\n",&NX);
    printf("Number of X intervals: %d.\n",NX);

    readString();
    fscanf(file_in,"%d\n",&NY);
    printf("Number of Y intervals: %d.\n",NY);

    NPX = NX+1;
    NPY = NY+1;

    readString();
    fscanf(file_in,"%d\n",&NTS);
    printf("Number of time steps: %d.\n",NTS);

    readString();
    fscanf(file_in,"%lf s\n",&dts);
    printf("Time step: %lf.\n",dts);

    readString();
    fscanf(file_in,"%d\n",&NI);
    printf("Maximum number of iterations per time step: %d.\n",NI);

    readString();
    fscanf(file_in,"%lf\n",&dt);
    dt = dts;
    printf("(USING DTS FOR NOW) Rate of convergence: %lf.\n",dt);

    readString();
    fscanf(file_in,"%lf\n",&urf);
    printf("Under Relaxation Factor (Pressure Correction): %lf.\n",urf);

    readString();
    fscanf(file_in,"%d\n",&control);
    printf("Time Step Control Check: %d.\n",control);

    int sw;
    double r, r_old, r_max;
    //double I_x, I_y;

    r_old = 1;

    readString();
    fscanf(file_in,"%d\n",&sw);
    printf("Number of sweeps: %d.\n",sw);

    readString();
    fscanf(file_in,"%lf\n",&r_max);
    printf("Maximum source term: %lf.\n",r_max);

    double XD, YD, U_in;

    readString();
    fscanf(file_in,"%lf\n",&XD);
    printf("X length: %lf.\n",XD);

    readString();
    fscanf(file_in,"%lf\n",&YD);
    printf("Y length: %lf.\n",YD);

    readString();
    fscanf(file_in,"%lf\n",&U_in);
    printf("Inlet speed: %lf.\n",U_in);

    double rho, vis;

    readString();
    fscanf(file_in,"%lf\n",&rho);
    printf("Density: %lf.\n",rho);

    readString();
    fscanf(file_in,"%lf\n",&vis);
    printf("Dynamic viscosity: %lf.\n",vis);

    double VAR1, VAR2, VAR3;

    readString();
    fscanf(file_in,"%lf\n",&VAR1);
    printf("Initial Temperature: %lf.\n",VAR1);

    readString();
    fscanf(file_in,"%lf\n",&VAR2);
    printf("Specific Heat: %lf.\n",VAR2);

    readString();
    fscanf(file_in,"%lf\n",&VAR3);
    printf("Thermal Conductivity: %lf.\n",VAR3);

    double cL, cS, cT;

    readString();
    fscanf(file_in,"%lf\n",&cL);
    printf("Characteristic Length: %lf.\n",cL);

    readString();
    fscanf(file_in,"%lf\n",&cS);
    printf("Characteristic Speed: %lf.\n",cS);

    readString();
    fscanf(file_in,"%lf\n",&cT);
    printf("Characteristic Time: %lf.\n",cT);

    fclose(file_in);

    double U_f, V_f;

    printf("\nPress any key to start calculations.");
    getchar();

    file = fopen("./Results/data.dat","w");
    file_R = fopen("./Results/mass_source.dat","w");

    //PID RELATED VARIABLES

    file_MDO = fopen("./Results/MDO.dat","w");
    //double PID_c = 0.01;

    //PID END

    int i, j, k, t, it;

    // GRID GENERATION

    //GRID(XD,YD,NPX,NPY,X,Y,XU,XV,YU,YV,DX,DY,DXU,DXV,DYU,DYV);
    //GRID_CAVITY(XD,YD,NPX,NPY,X,Y,XU,XV,YU,YV,DX,DY,DXU,DXV,DYU,DYV);
    GRID_BFS(XD,YD,NPX,NPY,X,Y,XU,XV,YU,YV,DX,DY,DXU,DXV,DYU,DYV);

    // PROBLEM INITIALIZATION

    INIT(U,NL,NC,0);
    INIT(V,NL,NC,0);
    INIT(U_NEW,NL,NC,0);
    INIT(V_NEW,NL,NC,0);
    INIT(P,NL,NC,101325);
    INIT(PC,NL,NC,0);
    INIT(R,NL,NC,0);

    INIT(RHO,NL,NC,rho);
    INIT(VIS,NL,NC,vis);
    INIT(VIS_T,NL,NC,0);

    INIT(T,NL,NC,VAR1);
    INIT(CP,NL,NC,VAR2);
    INIT(K,NL,NC,VAR3);

    INIT(I_x,NL,NC,0);
    INIT(I_y,NL,NC,0);

    fprintf(file,"VARIABLES='X','Y','U','V','P','T','Re','RHO',VIS','VIS_T','w','C','M'\n");

    //CALCULATION LOOP

    it = 0;

    for(t=0;t<=NTS;t++)
    {
        for(k=1;k<=NI;k++)
        {
            //BC(U,V,P,NPX,NPY,U_in,X,Y,t*dts,T,PID_c);
            //BC_CAVITY(U,V,P,NPX,NPY,U_in,X,Y,t*dts,T);
            BC_BFS(U,V,P,NPX,NPY,U_in,X,Y,t*dts,T);

            CALC_U(U,V,P,NPX,NPY,Y,YV,DX,DY,DXU,DXV,DYU,DYV,dt,RHO,VIS,VIS_T);
            CALC_V(U,V,P,NPX,NPY,X,XU,DX,DY,DXU,DXV,DYU,DYV,dt,RHO,VIS,VIS_T);
            CONTINUITY(U,V,RHO,R,NPX,NPY,DX,DY,DXU,DXV,DYU,DYV,dt);

            for(i=1;i<=sw;i++)
            {
                TDMA_X(PC,R,A,B,C,D,NPX,NPY,DX,DY,DXU,DXV,DYU,DYV,dt);
                TDMA_Y(PC,R,A,B,C,D,NPX,NPY,DX,DY,DXU,DXV,DYU,DYV,dt);
            }

            CALC_U_C(U,V,PC,RHO,NPX,NPY,DX,DY,DXU,DXV,DYU,DYV,dt);
            CALC_V_C(U,V,PC,RHO,NPX,NPY,DX,DY,DXU,DXV,DYU,DYV,dt);
            CALC_P(P,PC,NPX,NPY,urf);

            for(i=0;i<=NPX;i++)
            {
                for(j=0;j<=NPY;j++)
                {
                    U_f = (X[i]-XU[i])/(XU[i+1]-XU[i])*U[i+1][j]+(1-(X[i]-XU[i])/(XU[i+1]-XU[i]))*U[i][j];
                    V_f = (Y[j]-YV[j])/(YV[j+1]-YV[j])*V[i][j+1]+(1-(Y[j]-YV[j])/(YV[j+1]-YV[j]))*V[i][j];
                    Re[i][j] = RHO[i][j]*sqrt(U_f*U_f+V_f*V_f)*cL/VIS[i][j];
                }
            }

            PROPS(T,P,vis,VIS,VIS_T,RHO,CP,K,NPX,NPY,U,V,X,Y,DX,DY,DXU,DYU,DXV,DYV,I_x,I_y,Re,dts);

            Transform_Absolute_Element_Matrix(R,NPX,NPY);
            r = Maximum_Element_Matrix_NO_BOUNDARIES(R,NPX,NPY);
            it = it + 1;
            fprintf(file_R,"%d %lf \n",it,r);

            if ( t%control == 0 )
            {
                printf("\nTime Step: %d | Flow time: %lf | Iteration: %6.0d | Continuity Residual: %lf | Solution imbalance: %lf",t,t*dts,k,r,fabs(r-r_old)/fabs(r));
            }

            if ( (r < r_max && fabs(r-r_old)/fabs(r) <= 0.01) || fabs(r-r_old)/fabs(r) <= 0.001 )
            //if ( r < r_max || fabs(r-r_old)/fabs(r) <= 0.001 )
                break;

            //SOLUTION CONTROL
            if ( r > 100000 )
            {
                printf("\n\nFORCED EXIT DUE TO DIVERGENCE.");
                goto EXIT;
            }

            r_old = r;
        }

        //BC(U,V,P,NPX,NPY,U_in,X,Y,t*dts,T,PID_c);
        //BC_CAVITY(U,V,P,NPX,NPY,U_in,X,Y,t*dts,T);
        BC_BFS(U,V,P,NPX,NPY,U_in,X,Y,t*dts,T);

        if ( t%control == 0 )
        {
            //fprintf(file,"ZONE T='RESULTS' I=%d J=%d F=POINT C=BLACK SolutionTime=%d\n",NPX,NPY,t);
            fprintf(file,"ZONE T='RESULTS' I=%d J=%d F=POINT C=BLACK\n",NPX,NPY);

            for(j=0;j<NPY;j++)
            {
                for(i=0;i<NPX;i++)
                {
                    if ( i == 0 || i == NPX-1 || j == 0 || j == NPY-1 )
                    {
                        w = 0;
                        C1 = 0;
                    }
                    else
                    {
                        w = (V[i+1][j]-V[i][j])/DXV[i]-(U[i][j+1]-U[i][j])/DYU[j];
                        C1 = dts*(abs(U[i][j]/DXU[i])+abs(V[i][j]/DYV[j]));

                        if( C1 > 1.0 )
                        {
                            printf("\n\nALERT! Courant number exceeded 1.0!");
                            dts = dts*0.5;
                            printf("\nTime step reduction activated. New time step: %lfs",dts);
                            getchar();
                        }
                    }

                    U_f = (X[i]-XU[i])/(XU[i+1]-XU[i])*U[i+1][j]+(1-(X[i]-XU[i])/(XU[i+1]-XU[i]))*U[i][j];
                    V_f = (Y[j]-YV[j])/(YV[j+1]-YV[j])*V[i][j+1]+(1-(Y[j]-YV[j])/(YV[j+1]-YV[j]))*V[i][j];

                    //Re[i][j] = RHO[i][j]*sqrt(U_f*U_f+V_f*V_f)*cL/VIS[i][j];
                    Re[i][j] = RHO[i][j]*sqrt(U_f*U_f+V_f*V_f)*sqrt(DXU[i+1]*DXU[i+1]+DYV[j+1]*DYV[j+1])/VIS[i][j];

                    fprintf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",X[i],Y[j],U_f,V_f,P[i][j],T[i][j],Re[i][j],RHO[i][j],VIS[i][j],VIS_T[i][j],w,C1,R[i][j]);
                }
            }
        }

        //Next Time step

        //PID CALCULATION
        //PID_c = PID_c + PID(U,V,NPX,NPY,dts,U_in);
        //fprintf(file_MDO,"%lf %lf\n",t*dts,PID_c);

        PROPS(T,P,vis,VIS,VIS_T,RHO,CP,K,NPX,NPY,U,V,X,Y,DX,DY,DXU,DYU,DXV,DYV,I_x,I_y,Re,dts);
        CALC_U_t(U_NEW,P,U,V,NPX,NPY,Y,YV,DX,DY,DXU,DXV,DYU,DYV,dts,I_x,RHO,VIS,VIS_T);
        CALC_V_t(V_NEW,P,U,V,NPX,NPY,X,XU,DX,DY,DXU,DXV,DYU,DYV,dts,I_y,RHO,VIS,VIS_T);
        CALC_T_t(T_NEW,T,U,V,NPX,NPY,DX,DY,DXU,DXV,DYU,DYV,dts,RHO,VIS,CP,K);
    }

    //POS PROCESSING
    //LAW_WALL(vis,rho,NPX,NPY,U,V,X,Y,DX,DY,DXU,DYU,DXV,DYV);

    EXIT:

    fclose(file);
    fclose(file_MDO);

    for(i=0;i<NL;i++)
    {
        free(U[i]);
        free(U_NEW[i]);
        free(V[i]);
        free(V_NEW[i]);
        free(P[i]);
        free(PC[i]);
        free(R[i]);
        free(RHO[i]);
        free(VIS[i]);
        free(VIS_T[i]);
        free(T[i]);
        free(T_NEW[i]);
        free(CP[i]);
        free(K[i]);
        free(Re[i]);
        free(I_x[i]);
        free(I_y[i]);
    }

    free(U);
    free(U_NEW);
    free(V);
    free(V_NEW);
    free(P);
    free(PC);
    free(R);
    free(RHO);
    free(VIS);
    free(VIS_T);
    free(T);
    free(T_NEW);
    free(CP);
    free(K);
    free(Re);
    free(I_x);
    free(I_y);

    free(X);
    free(Y);
    free(XU);
    free(XV);
    free(YU);
    free(YV);
    free(DX);
    free(DY);
    free(DXU);
    free(DXV);
    free(DYU);
    free(DYV);

    free(A);
    free(B);
    free(C);
    free(D);

    printf("\n\nDone.\n");
    getchar();

    return 1;
}
