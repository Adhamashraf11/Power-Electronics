/*
 *  Name        : Power_Electronics.h
 *  Created on  : DEC 18, 2023
 *  Description : this file contains the application of bipolar selective harmonic elimination is used to calculate suitable angles
 *  Version     : v1.1
 *  Author      : Adham Ashraf Mohamed
 */

#include <stdio.h>
#include <math.h>

#define V_DC                 (      400       )
#define V_PEAK_FUNDAMENTAL   (      200       )
#define PI                   ( 3.14159265359  )
#define DEGREE_TO_RAD        (    PI/180.0    )
#define RAD_TO_DEGREE        (    180.0/PI    )
#define MAX_ITERATIONS       (       3        )
#define NUMBER_OF_ANGLES     (       6        )
#define NEEDED_Function      (NUMBER_OF_ANGLES)
#define JACOBIAN_SIZE_MATREX (NUMBER_OF_ANGLES)
#define CONST_X_V_DC         (((2*V_DC) / PI) )
#define DEBUGGER             (      0         )
    #if (DEBUGGER == 1)
        #define DEBUGG
    #else
        #undef DEBUGG
    #endif


void Get_Alphas             (int Number_Of_Angle     , double* InitialValue                               );
void Formalize_Function     (int Number_Of_Function  , double* InitialValue , double* Functions_Frame     );
void Get_Jacobian_Matrix    (int Size_Of_Matrex      , double* InitialValue , double (*Return_Jacobian)[JACOBIAN_SIZE_MATREX]);
void swap                   (double *Number_1        , double *Number_2);
void Inverse_Jacobian_Matrix(double (*Jacobian_Matrex_Frame)[JACOBIAN_SIZE_MATREX], double (*Return_Jacobian_Matrex_Inverse)[JACOBIAN_SIZE_MATREX] );
void Multiply_InvJacobian_Functions(double (*Inverse_Jacobian_Matrix)[JACOBIAN_SIZE_MATREX],double* Functions_Format,double *ReturnMultiply_Matrices);
void Subtract_Angles_MulMatrices(double *Old_Angle,double *Multiply_Matrices,double *Result);
void Print_Details (int Number_Of_Angle ,int Number_Of_iteration,double* Initial_Angles,double* New_Angles,double* Old_Angles,double (*Inverse_Jacobian_Matrix)[JACOBIAN_SIZE_MATREX],double* Functions_Format );
void Model_details (double* New_Angles,double* Old_Angles,double (*Inverse_Jacobian_Matrix)[JACOBIAN_SIZE_MATREX],double* Functions_Format);
void Subtract_Angles_MulMatrices_2(double *Old_Angle,double *Multiply_Matrices,double *Result_Degree);
int main(void)
{
    double InitialAngle[NUMBER_OF_ANGLES]={0}    ;  // array  to save alphas
    double Functions_Frame[NEEDED_Function]={0}  ;  // Model of functions_harmonics
    double Jacobian_Matrix[JACOBIAN_SIZE_MATREX][JACOBIAN_SIZE_MATREX];
    double Jacobian_Matrix_Inverse[JACOBIAN_SIZE_MATREX][JACOBIAN_SIZE_MATREX];
    double Multiply_Matrices[NUMBER_OF_ANGLES]={0};
    double Old_Angle[NUMBER_OF_ANGLES]={0}  ;
    double Result_Degree[NUMBER_OF_ANGLES]={0};
    int counter = 0 ;
    int ITERATIONS = 0;
     /* Get the initial Values of Angles */
    Get_Alphas(NUMBER_OF_ANGLES,InitialAngle);
        for(counter =0 ; counter < NUMBER_OF_ANGLES ; counter++)
        {
            Old_Angle[counter] =InitialAngle[counter] ;
        }

    for (ITERATIONS = 0 ; ITERATIONS < MAX_ITERATIONS ; ITERATIONS++)
    {
        if (ITERATIONS!=MAX_ITERATIONS)
        {
            /*Model of functions_harmonics to arry */
            Formalize_Function (NEEDED_Function,Old_Angle,Functions_Frame);

            /*Jacobian Matrix Frame*/
            Get_Jacobian_Matrix (JACOBIAN_SIZE_MATREX,Old_Angle,Jacobian_Matrix);

            /*Inversation Jacobian Matrix */
            Inverse_Jacobian_Matrix(Jacobian_Matrix,Jacobian_Matrix_Inverse);

            /*Multiply Inverse Jacobian Matrix with Functions format  */
            Multiply_InvJacobian_Functions(Jacobian_Matrix_Inverse,Functions_Frame,Multiply_Matrices);

            if (ITERATIONS==0)
            {
             /*Subtract Old Angles by Multiply Inverse Jacobian Matrix with Functions format  */
             Subtract_Angles_MulMatrices(Old_Angle,Multiply_Matrices,Result_Degree);
            }
            else
            {
              /*Subtract Old Angles by Multiply Inverse Jacobian Matrix with Functions format  */
             Subtract_Angles_MulMatrices_2(Old_Angle,Multiply_Matrices,Result_Degree);

            }
            }
        else
        {
            printf("MAX ITERATIONS !!!!!!!!!!!!!!!!!!!!!!!!!");
        }


        if ( Result_Degree[0]>100 )
        {
           printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
           printf("ERROR 404 !!!!!!!!!!!!!!!!!\n");
           printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

            break;
        }
        else
        {
            if (ITERATIONS == MAX_ITERATIONS-1)
            {
                 break ;
            }
            else
            {
                 for( counter =0 ; counter<NUMBER_OF_ANGLES; counter++)
                {
                Old_Angle[counter] = Result_Degree[counter] ;
                }
            }
        }

    }




#ifdef DEBUGG
#warning prints is closed
#else
     //  Print_Details(NUMBER_OF_ANGLES,MAX_ITERATIONS,InitialAngle, Result_Degree, Old_Angle, Jacobian_Matrix_Inverse , Functions_Frame );
        Model_details ( Result_Degree,Old_Angle,Jacobian_Matrix_Inverse,Functions_Frame);
#endif // DEBUGG

    return 0;
}




void swap(double *Number_1, double *Number_2)
{
    double temp = 0 ;

     temp     = *Number_1;
    *Number_1 = *Number_2;
    *Number_2 = temp;
}
void Get_Alphas(int Number_Of_Angle , double* InitialValue)
{
    int itera = 0 ;
    printf("Enter %d Angle Please\n", Number_Of_Angle);

    for (itera = 0; itera < Number_Of_Angle; itera++)
    {
        printf("Angle (%d) =", itera+1);
        scanf("%lf", &InitialValue[itera]);  // Corrected format specifier to %lf
    }

    #ifdef DEBUGG
        printf("debug Alphas Gauss (%d) = [ ",Number_Of_Angle);
        for (int debug1=0 ; debug1<Number_Of_Angle ; debug1++)
            {
                printf(" %f ",InitialValue[debug1]);
            }
            printf("]\n");
            printf("======================================================================================\n");
    #endif // DEBUGG

}
void Formalize_Function (int Number_Of_Function , double* InitialValue, double* Functions_Frame)
{
    int index = 0 ;
    int index_2 = 0;
    int n = 1 ; // Number of  Harmonics
    for (index_2 = 0 ; index_2 < Number_Of_Function ; index_2++)
    {
        if(index_2 == 0)
        {
           Functions_Frame[index_2 ] = ((CONST_X_V_DC/n)* (1 - (2 * cos(n * InitialValue[index+0] * DEGREE_TO_RAD))
                                                             + (2 * cos(n * InitialValue[index+1] * DEGREE_TO_RAD))
                                                             - (2 * cos(n * InitialValue[index+2] * DEGREE_TO_RAD))
                                                             + (2 * cos(n * InitialValue[index+3] * DEGREE_TO_RAD))
                                                             - (2 * cos(n * InitialValue[index+4] * DEGREE_TO_RAD))
                                                             + (2 * cos(n * InitialValue[index+5] * DEGREE_TO_RAD))))-V_PEAK_FUNDAMENTAL;

        }
        else
        {
            Functions_Frame[index_2]=(1 - (2 * cos(n * InitialValue[index+0] * DEGREE_TO_RAD))
                                        + (2 * cos(n * InitialValue[index+1] * DEGREE_TO_RAD))
                                        - (2 * cos(n * InitialValue[index+2] * DEGREE_TO_RAD))
                                        + (2 * cos(n * InitialValue[index+3] * DEGREE_TO_RAD))
                                        - (2 * cos(n * InitialValue[index+4] * DEGREE_TO_RAD))
                                        + (2 * cos(n * InitialValue[index+5] * DEGREE_TO_RAD)));
        }
         n=n+2;
    }

    #ifdef DEBUGG
        printf("debug function Matrix (%d) = [ ",Number_Of_Function);
        for (int debug1=0 ; debug1<NUMBER_OF_ANGLES ; debug1++)
            {
                printf(" %f ",Functions_Frame[debug1]);
            }
            printf("]\n");
            printf("-------------------------------------------------------------\n");
    #endif // DEBUGG


}
void Get_Jacobian_Matrix (int Size_Of_Matrex ,  double* InitialValue ,  double (*Return_Jacobian)[6])
{
    int row = 0 ;
    int col = 0 ;
    int n   = 1 ;
    for(row = 0 ; row < Size_Of_Matrex ; row++ )
    {
        for(col = 0 ; col < Size_Of_Matrex ; col++ )
        {

            if(row ==0)
            {
                if (col % 2 == 0)
                {       /*  EVEV COLUMS */
                    Return_Jacobian[row][col]=((CONST_X_V_DC/n)*(+2*n*sin(n*InitialValue[col]*DEGREE_TO_RAD)));
                }
                else
                {       /*  ODD COLUMS */
                    Return_Jacobian[row][col]=((CONST_X_V_DC/n)*(-2*n*sin(n*InitialValue[col]*DEGREE_TO_RAD)));
                }
            }
            else
            {
                if (col % 2 == 0)
                {       /*  EVEV COLUMS */
                    Return_Jacobian[row][col]=(+2*n*sin(n*InitialValue[col]*DEGREE_TO_RAD));
                }
                else
                {       /*  ODD COLUMS */
                    Return_Jacobian[row][col]=(-2*n*sin(n*InitialValue[col]*DEGREE_TO_RAD));
                }
            }
        }
        n+=2;
    }
    #ifdef DEBUGG
        printf("debug Jacobian Matrix (%d,%d) = \n",Size_Of_Matrex,Size_Of_Matrex);
        for (int debug1=0 ; debug1<Size_Of_Matrex ; debug1++)
            {
                for (int debug2=0 ; debug2<Size_Of_Matrex ; debug2++)
                {
                  printf(" %f ",Return_Jacobian[debug1][debug2]);
                }
             printf("\n");
            }
            printf("-------------------------------------------------------------\n");
    #endif // DEBUGG

}

void Inverse_Jacobian_Matrix(double (*Jacobian_Matrex_Frame)[JACOBIAN_SIZE_MATREX], double (*Return_Jacobian_Matrex_Inverse)[JACOBIAN_SIZE_MATREX])
 {
    // Initialize the inverse matrix as an identity matrix
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (i == j)
                Return_Jacobian_Matrex_Inverse[i][j] = 1.0;
            else
                Return_Jacobian_Matrex_Inverse[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 6; i++) {
        // Find pivot for column i
        int pivot = i;
        for (int j = i + 1; j < 6; j++) {
            if (Jacobian_Matrex_Frame[j][i] > Jacobian_Matrex_Frame[pivot][i])
                pivot = j;
        }

        // Swap rows i and pivot
        if (pivot != i) {
            for (int k = 0; k < 6; k++) {
                swap(&Jacobian_Matrex_Frame[i][k], &Jacobian_Matrex_Frame[pivot][k]);
                swap(&Return_Jacobian_Matrex_Inverse[i][k], &Return_Jacobian_Matrex_Inverse[pivot][k]);
            }
        }

        // Reduce rows
        double divisor = Jacobian_Matrex_Frame[i][i];
        for (int j = 0; j < 6; j++) {
            Jacobian_Matrex_Frame[i][j] /= divisor;
            Return_Jacobian_Matrex_Inverse[i][j] /= divisor;
        }

        for (int j = 0; j < 6; j++) {
            if (j != i) {
                double factor = Jacobian_Matrex_Frame[j][i];
                for (int k = 0; k < 6; k++) {
                    Jacobian_Matrex_Frame[j][k] -= factor * Jacobian_Matrex_Frame[i][k];
                    Return_Jacobian_Matrex_Inverse[j][k] -= factor * Return_Jacobian_Matrex_Inverse[i][k];
                }
            }
        }
    }
 #ifdef DEBUGG
        printf("debug Inverse Jacobian Matrix (%d,%d) = \n",JACOBIAN_SIZE_MATREX,JACOBIAN_SIZE_MATREX);
        for (int debug1=0 ; debug1<JACOBIAN_SIZE_MATREX ; debug1++)
            {
                for (int debug2=0 ; debug2<JACOBIAN_SIZE_MATREX; debug2++)
                {
                  printf(" %f ",Return_Jacobian_Matrex_Inverse[debug1][debug2]);
                }
             printf("\n");
            }
            printf("-------------------------------------------------------------\n");
    #endif // DEBUGG

}

void Multiply_InvJacobian_Functions(double (*Inverse_Jacobian_Matrix)[JACOBIAN_SIZE_MATREX],double* Functions_Format,double *ReturnMultiply_Matrices)
{

    int row = 0 ;
    int col = 0 ;

    for (row = 0; row < JACOBIAN_SIZE_MATREX; row++)
    {
        for (col = 0; col < JACOBIAN_SIZE_MATREX; col++)
        {
            ReturnMultiply_Matrices[row] =(ReturnMultiply_Matrices[row])+ ((Inverse_Jacobian_Matrix[row][col])*(Functions_Format[col]));
        }
    }
    #ifdef DEBUGG
        printf("debug Multiply Functions (%d) =[ ",JACOBIAN_SIZE_MATREX);
        for (int debug1=0 ; debug1<JACOBIAN_SIZE_MATREX ; debug1++)
            {
                printf(" %f ",ReturnMultiply_Matrices[debug1]);
            }
            printf("]\n");
            printf("-------------------------------------------------------------\n");
    #endif // DEBUGG

}

void Subtract_Angles_MulMatrices(double *Old_Angle,double *Multiply_Matrices,double *Result_Degree)
{

    int itera = 0 ;
    for( itera = 0 ; itera < NUMBER_OF_ANGLES ; itera++ )
    {
        Result_Degree[itera]=((Old_Angle[itera]*DEGREE_TO_RAD)- Multiply_Matrices[itera])*RAD_TO_DEGREE;
       // printf("New Angles =  %f \n",Result_Degree[itera]);
    }

    #ifdef DEBUGG
        printf("debug New Angle (%d) =[ ",NUMBER_OF_ANGLES);
        for (int debug1=0 ; debug1<NUMBER_OF_ANGLES ; debug1++)
            {
                printf(" %f ",Result_Degree[debug1]);
            }
            printf("]\n");
            printf("/*******************************************************************************************************************************/\n");
    #endif // DEBUGG

}

void Subtract_Angles_MulMatrices_2(double *Old_Angle,double *Multiply_Matrices,double *Result_Degree)
{

    int itera = 0 ;
    for( itera = 0 ; itera < NUMBER_OF_ANGLES ; itera++ )
    {
        Result_Degree[itera]=((Old_Angle[itera]*DEGREE_TO_RAD)- Multiply_Matrices[itera]*DEGREE_TO_RAD)*RAD_TO_DEGREE;
       // printf("New Angles =  %f \n",Result_Degree[itera]);
    }

    #ifdef DEBUGG
        printf("debug New Angle (%d) =[ ",NUMBER_OF_ANGLES);
        for (int debug1=0 ; debug1<NUMBER_OF_ANGLES ; debug1++)
            {
                printf(" %f ",Result_Degree[debug1]);
            }
            printf("]\n");
            printf("/*******************************************************************************************************************************/\n");
    #endif // DEBUGG

}


void Print_Details (int Number_Of_Angle ,int Number_Of_iteration,double* Initial_Angles,double* New_Angles,double* Old_Angles,double (*Inverse_Jacobian_Matrix)[JACOBIAN_SIZE_MATREX],double* Functions_Format )
{
    //int index = 0     ;
    int row   = 0     ;
    int col   = 0     ;
    printf("---------------------------------------------------------------------------------------------\n");
    printf("                Details List \n");
    printf("                ============ \n");

    printf("Number Of iteration = %d\n",Number_Of_iteration);
    printf("Number Of Angle     = %d\n",Number_Of_Angle);
    printf("V Peak Fundamental  = %d\n",V_PEAK_FUNDAMENTAL);

    printf("---------------------------------------------------------------------------------------------\n");
    printf("Initial  Angles = [ %f   , %f , %f , %f , %f , %f ]\n",Initial_Angles[0],Initial_Angles[1],Initial_Angles[2],Initial_Angles[3],Initial_Angles[4],Initial_Angles[5]);
    printf("Old  Angles     = [ %f   , %f , %f , %f , %f , %f ]\n",Old_Angles[0],Old_Angles[1],Old_Angles[2],Old_Angles[3],Old_Angles[4],Old_Angles[5]);
    printf("New  Angles     = [ %f   , %f , %f , %f , %f , %f ]\n",New_Angles[0],New_Angles[1],New_Angles[2],New_Angles[3],New_Angles[4],New_Angles[5]);
    printf("Functions_Frame = [ %f , %f , %f  , %f  , %f  , %f  ]\n",Functions_Format[0],Functions_Format[1],Functions_Format[2],Functions_Format[3],Functions_Format[4],Functions_Format[5]);
    printf("---------------------------------------------------------------------------------------------\n");

    for ( row=0 ; row<6 ; row++)
    {
        for (col=0 ; col<6 ; col++)
        {
            printf("Inverse Jacobian Matrix (%d,%d) = %f\n",row,col,Inverse_Jacobian_Matrix[row][col]);
        }
    }
    printf("---------------------------------------------------------------------------------------------\n");

}

void Model_details (double* New_Angles,double* Old_Angles,double (*Inverse_Jacobian_Matrix)[JACOBIAN_SIZE_MATREX],double* Functions_Format)
{

printf(" | %0.02f | ",New_Angles[0]);
printf("  |  %0.02f |",Old_Angles[0]);
printf("\t  | %0.04f %0.04f  %0.04f %0.04f  %0.04f %0.04f |",Inverse_Jacobian_Matrix[0][0],Inverse_Jacobian_Matrix[0][1],Inverse_Jacobian_Matrix[0][2],Inverse_Jacobian_Matrix[0][3],Inverse_Jacobian_Matrix[0][4],Inverse_Jacobian_Matrix[0][5]);
printf("\t| %0.04f |",Functions_Format[0]);
printf("\n");

printf(" | %0.02f | ",New_Angles[1]);
printf("  |  %0.02f |",Old_Angles[1]);
printf("\t  | %0.04f %0.04f %0.04f %0.04f  %0.04f %0.04f |",Inverse_Jacobian_Matrix[1][0],Inverse_Jacobian_Matrix[1][1],Inverse_Jacobian_Matrix[1][2],Inverse_Jacobian_Matrix[1][3],Inverse_Jacobian_Matrix[1][4],Inverse_Jacobian_Matrix[1][5]);
printf("\t| %0.04f   |",Functions_Format[1]);
printf("\n");


printf(" | %0.02f | =",New_Angles[2]);
printf(" |  %0.02f | -",Old_Angles[2]);
printf(" | %0.04f %0.04f %0.04f %0.04f %0.04f %0.04f |  * ",Inverse_Jacobian_Matrix[2][0],Inverse_Jacobian_Matrix[2][1],Inverse_Jacobian_Matrix[2][2],Inverse_Jacobian_Matrix[2][3],Inverse_Jacobian_Matrix[2][4],Inverse_Jacobian_Matrix[2][5]);
printf("\t|  %0.04f   |",Functions_Format[2]);
printf("\n");

printf(" | %0.02f |  ",New_Angles[3]);
printf(" |  %0.02f |",Old_Angles[3]);
printf("\t  | %0.04f %0.04f %0.04f %0.04f %0.04f %0.04f |",Inverse_Jacobian_Matrix[2][0],Inverse_Jacobian_Matrix[2][1],Inverse_Jacobian_Matrix[2][2],Inverse_Jacobian_Matrix[2][3],Inverse_Jacobian_Matrix[2][4],Inverse_Jacobian_Matrix[2][5]);
printf("\t|  %0.04f   |",Functions_Format[3]);
printf("\n");

printf(" | %0.02f |  ",New_Angles[4]);
printf(" |  %0.02f |",Old_Angles[4]);
printf("\t  | %0.04f %0.04f %0.04f %0.04f %0.04f %0.04f |",Inverse_Jacobian_Matrix[2][0],Inverse_Jacobian_Matrix[2][1],Inverse_Jacobian_Matrix[2][2],Inverse_Jacobian_Matrix[2][3],Inverse_Jacobian_Matrix[2][4],Inverse_Jacobian_Matrix[2][5]);
printf("\t|  %0.04f   |",Functions_Format[4]);
printf("\n");

printf(" | %0.02f |  ",New_Angles[5]);
printf(" |  %0.02f |",Old_Angles[5]);
printf("\t  | %0.04f %0.04f %0.04f %0.04f %0.04f %0.04f |",Inverse_Jacobian_Matrix[2][0],Inverse_Jacobian_Matrix[2][1],Inverse_Jacobian_Matrix[2][2],Inverse_Jacobian_Matrix[2][3],Inverse_Jacobian_Matrix[2][4],Inverse_Jacobian_Matrix[2][5]);
printf("\t|  %0.04f   |",Functions_Format[5]);
printf("\n");
printf("---------------------------------------------------------------------------------------------\n");


}







