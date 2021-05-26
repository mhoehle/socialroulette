// C++ code for solving the maximally diverse grouping problem.
//
// The algorithm is described in
// Xiangjing Lai and Jin-Kao Hao. "Iterated maxima search for the maximally
// diverse grouping problem. European Journal of Operational Research",
// 254(3), 780-800, 2016. https://doi.org/10.1016/j.ejor.2016.05.018
//
// A preprint is available from (Last query 2021-04-17):
// http://www.info.univ-angers.fr/%7Ehao/papers/LaiHaoEJOR2016.pdf
//
// Code by the authors is available from (Last query 2021-04-17):
// http://www.info.univ-angers.fr/pub/hao/mdgp.html
// and is copyright 2016 Xiangjing Lai and Jin-Kao Hao.
//
// The code was kindly donated for use in the package under a GPL v3 license
// The GPL v3 license terms are available from:
// https://www.gnu.org/licenses/gpl-3.0.txt or see the file LICENSE
//
// Minor modification of the code were performed in order to avoid R CMD check
// warnings (e.g. removing the definition of unused variables)

/******************************************************************************************/
/***************************    0. Head Files needed     **********************************/
/******************************************************************************************/
//#define WIN32_LEAN_AND_MEAN
#include <cstdlib>
#include <cstdlib>
#include "stdio.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <math.h>
#include <ctype.h>
#include <vector>
using namespace std;

#include <R.h>
#include <Rinternals.h>

/******************************************************************************************/
/********************    1. Data Structure and Global Variables    ************************/
/******************************************************************************************/
typedef struct Solution{
	int *p ;
	int *SizeG;
	double cost ;
}Solution ;

typedef struct Neighborhood{
	int  type ;
	int  v ;
    int  g ;
    int  x ;
    int  y ;
}Neighborhood;

char * File_Name;
char * Output_File_Name;
char * Solution_File;

int N, K;  // node number and group number
double f, f_best;
Solution CS, NS, GS, OS;
Neighborhood *Neighbors;
double total_time, starting_time, Time_limit;
int * p; // patition array for each vertex
int * bestp ; // patition array for each vertex
int * SizeG ;
double ** Delta_Matrix;  // incremental matrix
double ** Delta_Matrix1;
double ** Delta;
double ** Delta1;
double ** D;   // distance matrix between elements
double ** DT;
int * LB; // Lower bound of the number of elements in the group i
int * UB; // Upper bound of the number of elements in the group i

/******************************************************************************************/
/********************    2. Inputing Data and Assign Memeries   ***************************/
/******************************************************************************************/
void inputing()
{
    int i,j;
    int x1, x2;
    float d;
	ifstream FIC;
	FIC.open(File_Name);
	if ( FIC.fail() )
	{
		cout << "### Erreur open, File_Name " << File_Name << endl;
		exit(0);
	}
	FIC >> N >> K;
	char StrReading[100];
	FIC >> StrReading;
	if ( FIC.eof() )
	{
		cout << "### Error open, File_Name " << File_Name << endl;
		exit(0);
	}
	Rprintf("Read file...\n");
	if ( strcmp(StrReading, "ds" )==0 || strcmp(StrReading, "ss" )==0 )
	{
         LB = new int [K];
         UB = new int [K];
         for(i=0;i<K;i++) { FIC >> LB[i]; FIC >> UB[i];}
    }

    D = new double * [N];
    for(i=0;i<N;i++) D[i] = new double [N];
    for(i=0;i<N;i++)
      for(j=0;j<N;j++) D[i][j] = 0.0;

    DT = new double * [N];
    for(i=0;i<N;i++) DT[i] = new double [N];

    for(i=0;i<N;i++)
      for(j=0;j<N;j++) DT[i][j] = 0.0;

	while ( ! FIC.eof() )
	{
			FIC >> x1 >> x2 >> d;
			//cout << x1 <<"  "<< x2 <<"  "<<d<<" "<< endl;
			if ( x1<0 || x2<0 || x1>=N || x2 >=N )
			{
				cout << "### Error of node : x1="
					<< x1 << ", x2=" << x2 << endl;
				exit(0);
			}
			if(x1 != x2)
			{

                D[x2][x1] = d;
                D[x1][x2] = D[x2][x1] ;

                DT[x2][x1] = 2.0*d;
                DT[x1][x2] = DT[x2][x1] ;
			}

     }

     FIC.close();
 }

void AssignMemery()
{
    int i;

    p = new int [N];
    bestp = new int [N];
    SizeG = new int [K];

    Delta_Matrix = new double *[N];
       for(i=0;i<N;i++) Delta_Matrix[i] = new double [K];
    Delta_Matrix1 = new double *[N];
       for(i=0;i<N;i++) Delta_Matrix1[i] = new double [K];

    Delta = new double *[N];
       for(i=0;i<N;i++) Delta[i] = new double [K];

    Delta1 = new double *[N];
       for(i=0;i<N;i++) Delta1[i] = new double [K];

	CS.p = new int [N];
	NS.p = new int [N];
	GS.p = new int [N];
	OS.p = new int [N];

	CS.SizeG = new int [K];
	NS.SizeG = new int [K];
	GS.SizeG = new int [K];
	OS.SizeG = new int [K];

	Neighbors = new Neighborhood [N*(N-1)/2 + N*K ];
}

void ReleaseMemery()
{
     int i;

     delete [] p; p = NULL;
	 delete [] bestp; bestp = NULL;
	 delete [] SizeG; SizeG = NULL;

	 delete [] CS.p; CS.p = NULL;
	 delete [] CS.SizeG; CS.SizeG = NULL;
	 delete [] GS.p; GS.p = NULL;
	 delete [] GS.SizeG; GS.SizeG = NULL;
	 delete [] NS.p; NS.p = NULL;
	 delete [] NS.SizeG; NS.SizeG = NULL;
	 delete [] OS.p; OS.p = NULL;
	 delete [] OS.SizeG; OS.SizeG = NULL;

	 delete [] LB; LB = NULL;
	 delete [] UB; UB = NULL;
	 delete [] Neighbors; Neighbors = NULL;

	 for(i=0;i<N;i++)
	 {
	   delete [] Delta_Matrix[i]; Delta_Matrix[i] = NULL ;
       delete [] Delta_Matrix1[i]; Delta_Matrix1[i] = NULL ;
       delete [] Delta[i]; Delta[i] = NULL ;
       delete [] Delta1[i]; Delta1[i] = NULL ;
	   delete [] D[i]; D[i] = NULL;
       delete [] DT[i]; DT[i] = NULL;
	 }

}

/******************************************************************************************/
/*********************************    3. OutPuting Results   ******************************/
/******************************************************************************************/
int Proof(Solution &S)
{
    int i,j;
    double ff;
    int flag ;
    ff = 0.0;
    for( i = 0 ; i < N ; i++ )
	   for( j = i+1 ; j < N ; j++ )
	      {
	          if(S.p[i] == S.p[j])
               {
				  ff += D[i][j];
	           }
	      }
    S.cost = ff;
    for(i=0;i<K;i++) S.SizeG[i] = 0;
    for(i=0;i<N;i++) S.SizeG[S.p[i]]++;
    flag = 1;
    for(i=0; i < K; i++)
    if(S.SizeG[i] < LB[i]|| S.SizeG[i]> UB[i]) { flag = 0; break;}
   // Rprintf("  %d   ********** \n",flag);
    return flag;
}

void Outputing(Solution &S, char *filename)
{
  int i;int r;
	FILE *fp;
	char buff[80];
	r= rand()%1000;//random number in the range 0-1000

	if(Proof(S)==0) return;
    sprintf(buff,"%s",filename);
    fp=fopen(buff,"a+");
    fprintf(fp,"N = %d  G = %d  f = %lf\n", N , K, S.cost);
    for(i=0;i<K;i++)
    fprintf(fp,"%5.2d   %5.2d   %5.2d \n", LB[i], UB[i], S.SizeG[i]);
    //Rprintf("\n");
    for(i=0;i<N;i++)
    fprintf(fp,"%5.3d   %5.2d\n",i, S.p[i]);
	fclose(fp);
}

void Out_results(double best , double ave,  double worst, char *filename, char instance[])
{
	FILE *fp;
	char buff[80];
    sprintf(buff,"%s",filename);
    fp = fopen(buff,"a+");
    fprintf(fp,"%s   %lf   %lf   %lf\n", instance, best, ave, worst);
	fclose(fp);
}

/******************************************************************************************/
/****************   4. Constructive Heuristics for Initial Solution   *********************/
/******************************************************************************************/
void RandomInitiaSol(int p[],int SizeG[])
{

     int i;
     int p1;
     int count;
     int tot_number;
     int sum=0;
     int *Flag;
     int *SizeGroup;
     SizeGroup = new int [K];
     for(i=0;i<K;i++) SizeGroup[i] = 0;
     Flag = new int [N];
     for(i=0;i<N;i++) Flag[i] = 0;
     for(i=0; i < K; i++) sum += LB[i];
     tot_number=0;
     while(tot_number < sum)
     {
        p1 = rand()%N;
        if(Flag[p1]==0)
        {
           count=0;
           while(count<K)
           {
              if(SizeGroup[count] < LB[count])
              {
                p[p1] = count;
                Flag[p1] = 1;
                SizeGroup[count]++;
                tot_number++;
                break;
              }
              else count++;
           }
        }

     }

     tot_number=0;
     while(tot_number < N - sum)
     {
        p1 = rand()%N;
        if(Flag[p1]==0)
        {

           while(1)
           {
              count = rand()%K;
              if(SizeGroup[count] < UB[count])
              {
                p[p1] = count;
                Flag[p1] = 1;
                SizeGroup[count]++;
                tot_number++;
                break;
              }

           }
        }

     }

     for(i=0;i<K;i++)  SizeG[i] = SizeGroup[i];
     delete [] SizeGroup; SizeGroup= NULL;
     delete [] Flag; Flag = NULL;
}

void GreedyInitiaSol(int p[],int SizeGroup[])
{

     int i,j,g;
     int r, g_max;
     int cur_index;
     int tot_number;
     int sum = 0;
     int *Flag;
     double *SumG, MaxSumG;

     SumG = new double [K];
     for(g=0;g<K;g++) SizeGroup[g] = 0;
     Flag = new int [N];
     for(i=0;i<N;i++) Flag[i] = 0;

     //a. assign randomly K elements to K distinct groups
     for(g=0;g<K;g++)
     {
         while(1)
         {
             r = rand()%N;
             if(Flag[r]==0)
             {
                p[r] = g;
                Flag[r] = 1;
                SizeGroup[g]++;
                break;
             }
         }
     }

     //b. assign greedyly the elements to satisfy the lower bound constraints of groups
     for(g = 0; g < K; g++) sum += LB[g];
     tot_number = K;
     while(tot_number < sum)
     {

        cur_index = rand()%N;
        do
        {
           cur_index = (cur_index + 1)%N;
        }while(Flag[cur_index]);

        for(g=0;g<K;g++) SumG[g] = 0.0;

        for(j=0;j<N;j++)
           if(Flag[j]==1) SumG[p[j]] += D[cur_index][j];
        for(g=0;g<K;g++) SumG[g] /= SizeGroup[g];

        MaxSumG = -9999.0;
        for(g=0; g<K; g++)
         if(SizeGroup[g] < LB[g] && SumG[g] > MaxSumG)
         {
            MaxSumG = SumG[g];
            g_max = g;
         }

        p[cur_index] = g_max;
        Flag[cur_index] = 1;
        SizeGroup[g_max]++;

        tot_number++;
     }

     //c. assign the remaining the elements
     tot_number=0;
     while(tot_number < N - sum)
     {
        cur_index = rand()%N;
        do
        {
          cur_index = (cur_index + 1)%N;
        }while(Flag[cur_index]);

        for(g=0; g<K; g++) SumG[g] = 0.0;

        for(j=0;j<N;j++)
           if(Flag[j]==1) SumG[p[j]] += D[cur_index][j];
        for(g=0;g<K;g++) SumG[g] /= SizeGroup[g];

        MaxSumG = -9999.0;
        for(g=0; g<K; g++)
         if(SizeGroup[g] < UB[g] && SumG[g] > MaxSumG)
         {
            MaxSumG = SumG[g];
            g_max = g;
         }

        p[cur_index] = g_max;
        Flag[cur_index] = 1;
        SizeGroup[g_max]++;

        tot_number++;

    }

    delete [] Flag; Flag = NULL;
    delete [] SumG; SumG = NULL;
}

/******************************************************************************************/
/*******************************    5. Local Search Method   ******************************/
/******************************************************************************************/
void BuildNeighbors()
{
     int i,j,g;
     int count;
     //hoehle: not used - commented out: int SN = N*(N-1)/2 + N*K;
     count = 0;
     for(i=0;i<N;i++)
       for(g=0;g<K;g++)
       {
          Neighbors[count].type = 1;
          Neighbors[count].v = i ;
          Neighbors[count].g = g;
          count ++;
       }
     for(i=0;i<N;i++)
       for(j=i+1;j<N;j++)
       {
          Neighbors[count].type = 2;
          Neighbors[count].x = i;
          Neighbors[count].y = j;
          count ++;
       }
}

//2.1 Clear delta matrix
void Clear_Delta_Matrix( )
{
	int x, g ;
	f = 0.0;
	for( x = 0 ; x < N ; x++ )
		for( g = 0 ; g < K ; g++ )
		Delta_Matrix[ x ][ g ] = 0.0 ;
	//for( x = 0 ; x < N ; x++ )
	  //  for( g = 0 ; g < K ; g++ )
	  //  Delta[ x ][ g ] = 0.0 ;
	return ;
}

//2.2 Build delta matrix
void Build_Delta_Matrix()
{
	int i,j;
	Clear_Delta_Matrix( );

	for(i = 0; i < N ; i++ )
	   for( j = 0 ; j < N ; j++ )
	     {
		     Delta_Matrix[ i ][ p[j] ]  +=  D[i][j];
	     }

	for(i = 0; i < N ; i++ )
	   for( j = 0 ; j < K ; j++ )
       Delta[i][j] = Delta_Matrix[ i ][ j ] - Delta_Matrix[ i ][ p[i] ];

    f = 0.0;
    for( i = 0 ; i < N ; i++ )
      f += Delta_Matrix[i][p[i]];
    f = f/2;

	//printf("\n f= %lf ********** \n", f);
    return ;
}

//2.2 Update one move delta matrix
void One_Move_Update_Delta_Matrix(int i, int g0, int g1)
{
    int x,j,k;

    for(j=0;j<N;j++)
     {
           if(j!=i)
           {
             Delta_Matrix[ j ][ g0 ] -= D[i][j];
             Delta_Matrix[ j ][ g1 ] += D[i][j];
           }
     }

    for(x=0; x<N; x++)
     {
         if(x!=i)
           {
               Delta[x][g0] =   Delta_Matrix[ x ][ g0 ] - Delta_Matrix[ x ][ p[ x ] ];
               Delta[x][g1] =   Delta_Matrix[ x ][ g1 ] - Delta_Matrix[ x ][ p[ x ] ];

               if(p[x]==g0)
               {
                 for(k=0;k<K;k++)
                  {
                     Delta[x][k] =  Delta_Matrix[ x ][ k ] - Delta_Matrix[ x ][ g0 ];
                  }
               }

               if(p[x]==g1)
               {
                 for(k=0;k<K;k++)
                 {
                    Delta[x][k]=  Delta_Matrix[ x ][ k ] - Delta_Matrix[ x ][ g1 ];
                 }
               }

           }
     }
    x = i;
    for(k=0; k<K; k++) Delta[x][k] = Delta_Matrix[ x ][ k ] - Delta_Matrix[ x ][ g1 ];

	return ;
}
void One_Move_Update_Delta_Matrix1(int i, int g0, int g1)
{
    int j;

    for(j=0;j<N;j++)
     {
           if(j!=i)
           {
             Delta_Matrix[ j ][ g0 ] -= D[i][j];
             Delta_Matrix[ j ][ g1 ] += D[i][j];
           }
     }
	return ;
}

void RandLS(int partition[], int SizeGroup[], double *cost)
{
    int i, v, g, x, y;
    int old_g, old_g1, swap ;
    double delt;
    int Flag;

    for(i=0;i<N;i++) p[i] = partition[i];
	Build_Delta_Matrix();
	f_best = f ;

    delt = -99999.0;
    do
	{
        Flag = 0;

        for(v=0;v<N;v++)
        for(g=0;g<K;g++)
          if( ( p[v] != g ) && ( SizeGroup[ p[v] ] > LB[ p[v] ] ) && (SizeGroup[g] < UB[g])  )
             {
               delt = Delta_Matrix[ v ][ g ] - Delta_Matrix[ v ][ p[ v ] ];
               if(delt > 0)
               {
                  old_g = p[v] ;
				  One_Move_Update_Delta_Matrix1(v, old_g, g);
				  SizeGroup[old_g] = SizeGroup[old_g] - 1;
                  SizeGroup[g] = SizeGroup[g] + 1;
				  p[v] = g;
			  	  f += delt;
			  	  Flag = 1;
               }

             }

         for(x=0;x<N;x++)
            for(y=x+1;y<N;y++)
             if(p[x] != p[y])
             {
               delt = (Delta_Matrix[ x ][ p[y] ] - Delta_Matrix[ x ][ p[x] ]) + (Delta_Matrix[ y ][ p[x] ] - Delta_Matrix[ y ][ p[y] ]) - DT[x][y];
               if(delt > 0)
               {
                old_g  = p[ x ];
                old_g1 = p[ y ];
				One_Move_Update_Delta_Matrix1( x, old_g, old_g1);
				One_Move_Update_Delta_Matrix1( y, old_g1, old_g);

				swap   = p[ x ];
				p[ x ] = p[ y ];
				p[ y ] = swap;

				f += delt;
				Flag = 1;
               }

             }
     }while(Flag == 1);

    for(i=0;i<N;i++)  partition[i] = p[i];
	*cost = f;
//	printf("f=%lf\n",f);
}

/******************************************************************************************/
/**********************     6. Maxima Search  Procedure     *******************************/
/******************************************************************************************/
void MinimaSearch(int partition[], int SizeGroup[], double *cost)
{
    int i, j, v, g, x, y, L ;  int Q, Flag; int tt=0;


    int NumberNeighbors, old_g, old_g1, swap ;

    //hoehle: not used: int starting_index;
    int cur_index;
    int Counter, CC;
    double delt, delt_max; //hoehle: not used: delt1,

    NumberNeighbors = N*(N-1)/2 + N*K ;

    for(i=0;i<N;i++) p[i] = partition[i];
    for(j=0;j<K;j++) SizeG[j] = SizeGroup[j];

	Build_Delta_Matrix();
	f_best = f;
	Counter = 0;

    do
    {

      do
	  {
        Flag = 0;

        for(v=0;v<N;v++)
          for(g=0;g<K;g++)
           if( ( p[v] != g ) && ( SizeG[ p[v] ] > LB[ p[v] ] ) && (SizeG[g] < UB[g])  )
             {
               delt = Delta_Matrix[ v ][ g ] - Delta_Matrix[ v ][ p[ v ] ];
              // delt = Delta[v][g];
              // if(delt != delt1) printf("! errer delta one move ! \n");
              // else printf("equal! ");
              // if( Delta[v][g] > 0.0 )

               if(delt > 0)
               {
                  old_g = p[v] ;
				  One_Move_Update_Delta_Matrix(v, old_g, g);
				  SizeG[old_g] = SizeG[old_g] - 1;
                  SizeG[g] = SizeG[g] + 1;
				  p[v] = g;
			  	  f += delt;
			  	 // f += Delta[v][g];
			  	  Flag = 1;
               }
             }

         for(x=0;x<N;x++)
           for(y=x+1;y<N;y++)
             if(p[x] != p[y])
             {

             //  delt1 = (Delta_Matrix[ x ][ p[y] ] - Delta_Matrix[ x ][ p[x] ]) + (Delta_Matrix[ y ][ p[x] ] - Delta_Matrix[ y ][ p[y] ]) - DT[x][y];
               delt = Delta[ x ][ p[y] ] + Delta[ y ][ p[x] ] - DT[x][y];
             //  if(delt != delt1) printf("errer delta swap move!\n");

               if(delt > 0)
               {
                old_g  = p[ x ];
                old_g1 = p[ y ];
				One_Move_Update_Delta_Matrix( x, old_g, old_g1);

				swap   = p[ x ];
				p[ x ] = p[ y ];
				One_Move_Update_Delta_Matrix( y, old_g1, old_g);

				p[ y ] = swap;

				f += delt;
				Flag = 1;
               }
             }


       }while(Flag == 1);

       if(f > f_best)
       {
            f_best = f;
            *cost = f_best;
            for(j=0; j<N; j++) partition[j] = p[j];
            for(j=0; j<K; j++) SizeGroup[j] = SizeG[j];

            for(i=0;i<N;i++)
               for(j=0;j<K;j++) Delta_Matrix1[i][j] = Delta_Matrix[i][j];
            for(i=0;i<N;i++)
               for(j=0;j<K;j++) Delta1[i][j] = Delta[i][j];
            Counter = 0;
       }
       else Counter++;

       f = *cost;
       for(j=0; j<N; j++) p[j] = partition[j];
       for(j=0; j<K; j++) SizeG[j] = SizeGroup[j];
       for(i=0;i<N;i++)
          for(j=0;j<K;j++) Delta_Matrix[i][j] = Delta_Matrix1[i][j];
       for(i=0;i<N;i++)
          for(j=0;j<K;j++) Delta[i][j] = Delta1[i][j];

       for(L = 0; L < 3; L ++)
       {
         delt_max = -9999999;
         CC = 0;
         do
         {

           cur_index = rand()%NumberNeighbors; // Please ensure that MAX_RAND is large enough !!!!

           if(Neighbors[cur_index].type == 1 )
           {
             v =  Neighbors[cur_index].v;
             g =  Neighbors[cur_index].g;
             if( ( p[v] != g ) && ( SizeG[ p[v] ] > LB[ p[v] ] ) && (SizeG[g] < UB[g])  )
             {
                  delt = Delta_Matrix[ v ][ g ] - Delta_Matrix[ v ][ p[ v ] ];
                  if(delt > delt_max)  { Q = cur_index; delt_max = delt;  }
                  CC++;
             }
           }

           if(Neighbors[cur_index].type == 2 )
           {
             x =  Neighbors[cur_index].x;
             y =  Neighbors[cur_index].y;
             if(p[x] != p[y])
             {
               // delt = (Delta_Matrix[ x ][ p[y] ] - Delta_Matrix[ x ][ p[x] ]) + (Delta_Matrix[ y ][ p[x] ] - Delta_Matrix[ y ][ p[y] ]) - 2*D[x][y];
                delt = Delta[ x ][ p[y] ] + Delta[ y ][ p[x] ] - DT[x][y];
              //  if(delt != delt) printf("errer !\n");
                if(delt > delt_max)  { Q = cur_index; delt_max = delt; }
                CC++;
             }
           }

            //printf("cc=%d delt_max =%lf \n",CC, delt_max);
          }while(CC < N);

          if(Neighbors[Q].type == 1 )
          {
             v =  Neighbors[Q].v;
             g =  Neighbors[Q].g;
             delt = Delta_Matrix[ v ][ g ] - Delta_Matrix[ v ][ p[ v ] ];
             old_g = p[v] ;
		     One_Move_Update_Delta_Matrix(v, old_g, g);
		     SizeG[old_g] = SizeG[old_g] - 1;
             SizeG[g] = SizeG[g] + 1;
		     p[v] = g;
	  	     f += delt;
          }
          else if(Neighbors[Q].type == 2 )
           {
             x =  Neighbors[Q].x;
             y =  Neighbors[Q].y;
             delt = (Delta_Matrix[ x ][ p[y] ] - Delta_Matrix[ x ][ p[x] ]) + (Delta_Matrix[ y ][ p[x] ] - Delta_Matrix[ y ][ p[y] ]) - 2*D[x][y];

             old_g  = p[ x ];
             old_g1 = p[ y ];
	         One_Move_Update_Delta_Matrix( x, old_g, old_g1);

	         swap  = p[ x ];
	         p[ x ] = p[ y ];
	         One_Move_Update_Delta_Matrix( y, old_g1, old_g);

	         p[ y ] = swap;

	         f += delt;
           }

       }
      tt++;

    }while(Counter < 3);
   // printf("tt=%d\n",tt);
}

/******************************************************************************************/
/**********************  7. Strong Perturbation Operators     *****************************/
/******************************************************************************************/
void Shake1(int partition[], int k_max)
{
   int r1,r2;
   int swap;
   int iter;

   iter = k_max;
   while( iter > 0 )
   {
       r1 = rand()%N;
       while(1)
       {
          r2 = rand()%N;
          if( partition[r1] != partition[r2]) break;
       }
       swap = partition[r1];
       partition[r1] = partition[r2];
       partition[r2] = swap;

       iter--;
   }
}

void Shake2(int L, int partition[], int SizeGroup[])
{
    int i, v, g, x, y;
    int NumberNeighbors, old_g, old_g1, swap ;
    //hoehle: int iter=0;
    int cur_index; //hoehle:, starting_index;
    double delt;
    //hoehle - these are already global vars?
    //double total_time=0.0, starting_time=0.0;
    total_time=0.0, starting_time=0.0;
    int theta, count = 0;

    theta =  L;

    NumberNeighbors = N*(N-1)/2 + N*K ;
    for(i=0;i<N;i++) p[i] = partition[i];

    do
	{

        cur_index = rand()% NumberNeighbors;

        if(Neighbors[cur_index].type == 1 )
        {
             v =  Neighbors[cur_index].v;
             g =  Neighbors[cur_index].g;
             if( ( p[v] != g ) && ( SizeGroup[ p[v] ] > LB[ p[v] ] ) && (SizeGroup[g] < UB[g])  )
             {
               delt = Delta_Matrix[ v ][ g ] - Delta_Matrix[ v ][ p[ v ] ];

               old_g = p[v] ;
				  //One_Move_Update_Delta_Matrix(v, old_g, g);
               SizeGroup[old_g] = SizeGroup[old_g] - 1;
               SizeGroup[g] = SizeGroup[g] + 1;
               p[v] = g;
			  	  //f += delt;
	  	       count ++;
			  	 // printf("\n One_Move :  %8d       %5.6lf       %5.6lf       %5.3lf s", count, f, f_best, total_time);

             }

        }

        else if(Neighbors[cur_index].type == 2 )
        {
             x =  Neighbors[cur_index].x;
             y =  Neighbors[cur_index].y;
             if(p[x] != p[y])
             {
                delt = (Delta_Matrix[ x ][ p[y] ] - Delta_Matrix[ x ][ p[x] ]) + (Delta_Matrix[ y ][ p[x] ] - Delta_Matrix[ y ][ p[y] ]) - 2*D[x][y];

                old_g  = p[ x ];
                old_g1 = p[ y ];
				//One_Move_Update_Delta_Matrix( x, old_g, old_g1);
				//One_Move_Update_Delta_Matrix( y, old_g1, old_g);

				swap   = p[ x ];
				p[ x ] = p[ y ];
				p[ y ] = swap;

				//f += delt;
				count ++;
			 //	printf("\n One_Move :  %8d       %5.6lf       %5.6lf       %5.3lf s", count, f, f_best, total_time);
             }
        }
      //iter++;
   	  //if(f > f_best) f_best = f;
      //total_time = ( clock() - starting_time )/CLOCKS_PER_SEC ;

    }while(count < theta);

    for(i=0;i<N;i++)  partition[i] = p[i];
 }

/******************************************************************************************/
/***********************************    8. Initial solution    ****************************/
/******************************************************************************************/
void InitialSol(Solution &S)
{
     int i,j;
     int counter=0;
     GreedyInitiaSol(S.p, S.SizeG);
     RandLS(S.p, S.SizeG, &S.cost);
     counter++;
     while(counter < 10)
     {
        GreedyInitiaSol(NS.p, NS.SizeG);
        RandLS(NS.p, NS.SizeG, &NS.cost);
        if(NS.cost > S.cost)
        {
           for(i=0;i<N;i++) S.p[i] = NS.p[i];
           for(j=0;j<K;j++) S.SizeG[j] = NS.SizeG[j];
           S.cost = NS.cost;
        }
       counter++;
     }
}

/******************************************************************************************/
/*****************      9. Iterated Maxima Search Algorithm *******************************/
/******************************************************************************************/

void IMS()
{
    int i,j;
    int L;

    starting_time = clock();
    GS.cost = -9999999;
    InitialSol(CS);
    for(i=0;i<N;i++) GS.p[i] = CS.p[i];
    for(j=0;j<K;j++) GS.SizeG[j] = CS.SizeG[j];
    GS.cost = CS.cost;
    if(N<=400) L = int (1.0*N/K);
    else L = int (1.5*N/K);

    while(1.0*(clock()- starting_time)/CLOCKS_PER_SEC  <  Time_limit)
    {
      Shake2(L, CS.p, CS.SizeG);
      MinimaSearch(CS.p, CS.SizeG, &CS.cost);

      if(CS.cost > GS.cost)
      {
        for(i=0;i<N;i++) GS.p[i] = CS.p[i];
        for(j=0;j<K;j++) GS.SizeG[j] = CS.SizeG[j];
        GS.cost = CS.cost;
        //total_time = ( clock() - starting_time )/CLOCKS_PER_SEC ;
        //printf("f = %lf  time=%lf s \n", GS.cost,  total_time);
      }

    }
}

/******************************************************************************************/
/********************    10. Main Function for Multiple Runs     **************************/
/******************************************************************************************/

// int main(int argc, char *argv[])
// {
//     int i,j;
//     int i1,j1;
//     int seed;
//     const int  Times = 10; //hoehle: used to be 20
//     double F[Times];
//     double F_best= -99999999, F_worst = 999999999, F_ave = 0.0;
//     seed = time(NULL) % 1000000 ;
//      srand( seed );
//
//     /*
//     File_Name = "MDG-a_21.txt";
//     Output_File_Name = "ss.txt";
//     Solution_File = "MDG-a_21.sol";
//     Time_limit = 1200;
//     */
//
//     File_Name = argv[1];
//     Solution_File = argv[2];
//     Output_File_Name = argv[3];
//     Time_limit = atoi(argv[4]);
//     //hoehle: Output_File_Name = "new.txt";
//
//     Rprintf("Inputting...\n");
//     inputing();
//     Rprintf("Assigning memory...\n");
//     AssignMemery();
//
//     if(N==120) Time_limit = 3;
//     else if(N==240)Time_limit = 20;
//     else if(N==480)Time_limit = 120;
//     else if(N==960)Time_limit = 600;
//     else if(N==2000)Time_limit= 1200;
//     else if(N==3000)Time_limit = 3000;
//
//     printf("Building neighbours...\n");
//     BuildNeighbors();
//     printf("Running (this might take a while)...\n");
//     OS.cost = -99999.0;
//     for(j=0;j<Times;j++) F[j] = 0.0;
//     for(i=0; i < Times; i++)
//     {
//       IMS();
//       if(Proof(GS))
//       {
//         F[i] = GS.cost;
//         if(F[i]> OS.cost)
//         {
//             for(i1=0;i1<N;i1++) OS.p[i1] = GS.p[i1];
//             for(j1=0;j1<K;j1++) OS.SizeG[j1] = GS.SizeG[j1];
//             OS.cost = GS.cost;
//         }
//       }
//       printf("%lf \n", F[i]);
//     }
//     for(i=0;i<Times;i++)
//     {
//        if(F[i] > F_best )  F_best = F[i];
//        if(F[i] < F_worst)  F_worst = F[i];
//        F_ave += F[i];
//     }
//     F_ave /=  Times;
//     Out_results(F_best , F_ave, F_worst, Output_File_Name, File_Name);
//     Outputing(OS, Solution_File);
//     ReleaseMemery();
//     return 0;
// }

// Function for calling the solver from R
extern "C" {

  void mdgp(char **File_Name_R, char **Output_File_Name_R, char **Solution_File_Name_R, double *Time_limit_R) {
    int i,j;
    int i1,j1;
    int seed;
    const int  Times = 10; //hoehle: used to be 20
    double F[Times];
    double F_best= -99999999, F_worst = 999999999, F_ave = 0.0;

    //Seed the random number generator
    seed = time(NULL) % 1000000 ;
    srand( seed );
    //Setup R random number generator
    GetRNGstate();

    //Assign variables
    File_Name = *File_Name_R;
    Output_File_Name = *Output_File_Name_R; //"output.txt";
    Solution_File = *Solution_File_Name_R;
    Time_limit = *Time_limit_R;


    //Output variables
    Rprintf("Input file: %s\n", File_Name);
    Rprintf("Output file: %s\n", Output_File_Name);
    Rprintf("Solution file: %s\n", Solution_File);
    Rprintf("Time_limit: %f\n", Time_limit);


    Rprintf("Inputting...\n");
    inputing();
    Rprintf("Assigning memory...\n");
    AssignMemery();

    //Adjust time limit for special cases
    // if(N==120) Time_limit = 3;
    // else if(N==240)Time_limit = 20;
    // else if(N==480)Time_limit = 120;
    // else if(N==960)Time_limit = 600;
    // else if(N==2000)Time_limit= 1200;
    // else if(N==3000)Time_limit = 3000;

    Rprintf("Building neighbours...\n");
    BuildNeighbors();
    Rprintf("Running...\n");
    OS.cost = -99999.0;
    for(j=0;j<Times;j++) F[j] = 0.0;
    for(i=0; i < Times; i++)
    {
      IMS();
      if(Proof(GS))
      {
        F[i] = GS.cost;
        if(F[i]> OS.cost)
        {
          for(i1=0;i1<N;i1++) OS.p[i1] = GS.p[i1];
          for(j1=0;j1<K;j1++) OS.SizeG[j1] = GS.SizeG[j1];
          OS.cost = GS.cost;
        }
      }
      Rprintf("%lf \n", F[i]);
    }
    for(i=0;i<Times;i++)
    {
      if(F[i] > F_best )  F_best = F[i];
      if(F[i] < F_worst)  F_worst = F[i];
      F_ave += F[i];
    }
    F_ave /=  Times;

    Rprintf("Outputting...\n");
    Out_results(F_best , F_ave, F_worst, Output_File_Name, File_Name);
    Rprintf("Outputting solution file...\n");
    Outputing(OS, Solution_File);
    Rprintf("Done...cleaning memory and whatnot...\n");
    ReleaseMemery();

    //Release RNG - https://cran.r-project.org/doc/manuals/R-exts.html#Random-numbers
    PutRNGstate();

    //Done
  } /* end of mdgp */
} /* end of extern "C" */
