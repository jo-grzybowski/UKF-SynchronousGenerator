#ifndef HEADER3_H_
#define HEADER3_H_
#include<math.h>


const int n = 7; // Numero de estados do sistema + parametros
const int n_me = 4; // Numero de estados medidos
const int n_orig = 7; // Número de estados do sistema
const double Ts = 0.0001; // Tempo de amostragem
const double h = 0.00001; // Passo de integração
const double tf = 0.1;
//int N = tf/h+1;
const int N = 1001;
//int N_Ts = tf/Ts+1;
const int N_Ts = 1001;
const double Poles = 4.0; // número de polos do gerador
const double pi = 3.1415926535897;


// RANDOM VECTOR GEN 
void rdn_vec(double elements[], int nrows){
	int i;
	srand( (unsigned)time(NULL) );
	for (i = 0; i<nrows; i++) {
  		// generate random index
  		elements[i] = ((rand() % 100)/100.)*pow(-1,rand() % 9);
	};
}

// MULTIPLICAÇÃO MATRIZ-ESCALAR **************************************************************************************
void me_mult(double A[][7],int nrows,int ncols, double c){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			A[i][j] = c*A[i][j];
		}
	}
}

// MULTIPLICAÇÃO MATRIZ-ESCALAR **************************************************************************************
void me_mult_me(double A[][4],int nrows,int ncols, double c){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			A[i][j] = c*A[i][j];
		}
	}
}

// PRODUTO ESCALAR
double dot_p(double x[], double y[], int n)
{
    double res = 0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

// ADICAO DE VETORES
void add(double a[], double b[], int rows, double c[]){
	int i;
		for (i = 0;i<rows;i++){
			c[i] = a[i] + b[i];
	}
}

//SUBTRACAO DE VETORES
void subtr(double a[], double b[], int rows, double c[]){
	int i;
		for (i = 0;i<rows;i++){
			c[i] = a[i] - b[i];
	}
}

// ATRIBUIÇÃO
void atr(double a[], int rows, double c[]){
	int i;
		for (i = 0;i<rows;i++){
			c[i] = a[i];
	}
}

// MULTIPLICAÇÃO POR ESCALAR
void sc_vec(double a, double b[], int rows, double c[]){
	int i;
		for (i = 0;i<rows;i++){
			c[i] = a*b[i];
	}
}

// ADIÇÃO DE MATRIZES **************************************************************************************
void mm_add(double A[][7], double B[][7], int nrows, int ncols, double C[][7]){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

// ADIÇÃO DE MATRIZES **************************************************************************************
void mm_subtr(double A[][7], double B[][7], int nrows, int ncols, double C[][7]){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			C[i][j] = A[i][j] - B[i][j];
		}
	}
}

// ADIÇÃO DE MATRIZES **************************************************************************************
void mm_add_me(double A[][4], double B[][4], int nrows, int ncols, double C[][4]){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

// ADIÇÃO DE MATRIZES 2n+1
void mm_add_2np1(double A[][15], double B[][15], int nrows, int ncols, double C[][15]){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

// MATRIZ IDENTIDADE [n]
void eye(double A[][7], int ncols){
	int i,j;
	for (i = 0;i<ncols;i++){
		for (j = 0;j<ncols;j++){
			if (i == j){
				A[i][j] = 1.0;
			}
			if (i!=j){
				A[i][j] = 0.0;
			}
			
		}
	}
}

// MATRIZ IDENTIDADE [n_me]
void eye_me(double A[][4], int ncols){
	int i,j;
	for (i = 0;i<ncols;i++){
		for (j = 0;j<ncols;j++){
			if (i == j){
				A[i][j] = 1.0;
			}
			if (i!=j){
				A[i][j] = 0.0;
			}
			
		}
	}
}

// MATRIZ DE ZEROS [][n]
void zeros_matrix(double A[][7],int order){
	int i,j;
	for (i=0;i<order;i++){
		for (j=0;j<n;j++){
			A[i][j] = 0.0;
		}
	}
}

// MATRIZ DE ZEROS [][2n+1]
void zeros_matrix2np1(double A[][15],int order){
	int i,j;
	for (i=0;i<order;i++){
		for (j=0;j<2*order+1;j++){
			A[i][j] = 0.0;
		}
	}
}

// VETOR DE ZEROS
void zeros_vec(double a[], int size){
	int i;
	for (i=0;i<size;i++){
			a[i] = 0.0;
	}
}

// MULTIPLICAÇÃO MATRIZ-VETOR
void mv_mult(double A[][7], double x[], int rows, int cols, double b[]){
	int i,j;
    for (i = 0; i < rows; i++)
    {
    	b[i] = 0.0;
    	for (j = 0;j<cols;j++){
			b[i] += (A[i][j] * x[j]);
    }
    }
}

// MULTIPLICAÇÃO MATRIZ-VETOR TAMANHO ORIGINAL
void mv_mult_orig(double A[][7], double x[], int rows, int cols, double b[]){
	int i,j;
    for (i = 0; i < rows; i++)
    {
    	b[i] = 0.0;
    	for (j = 0;j<cols;j++){
			b[i] += (A[i][j] * x[j]);
    }
    }
}

// MULTIPLICAÇÃO MATRIZ-VETOR ME
void mv_mult_me(double A[][4], double x[], int rows, int cols, double b[]){
	int i,j;
    for (i = 0; i < rows; i++)
    {
    	b[i] = 0.0;
    	for (j = 0;j<cols;j++){
			b[i] += (A[i][j] * x[j]);
    }
    }
}

// MULTIPLICAÇÃO MATRIZ-VETOR 2n+1
void mv_mult_2np1(double A[][15], double x[], int rows, int cols, double b[]){
	int i,j;
    for (i = 0; i < rows; i++)
    {
    	b[i] = 0.0;
    	for (j = 0;j<cols;j++){
			b[i] += (A[i][j] * x[j]);
    }
    }
}

// ATRIB COLUMN
void atr_col(double A[][7], int ncols, int nrows, double vec[], int opt){
	// opt = 1 matriz -> vetor
	// opt = 2 vetor -> matriz
	int i;
	if (opt == 1){
		for (i = 0;i<nrows;i++){
			vec[i] = A[i][ncols];
		}
	}
	if (opt == 2){
		for (i = 0;i<nrows;i++){
			A[i][ncols] = vec[i];
		}
	}
}

// ATRIB COLUMN
void atr_col_me(double A[][4], int ncols, int nrows, double vec[], int opt){
	// opt = 1 matriz -> vetor
	// opt = 2 vetor -> matriz
	int i;
	if (opt == 1){
		for (i = 0;i<nrows;i++){
			vec[i] = A[i][ncols];
		}
	}
	if (opt == 2){
		for (i = 0;i<nrows;i++){
			A[i][ncols] = vec[i];
		}
	}
}

// ATRIB COLUMN N
void atr_col_N(double A[][1001], int ncols, int nrows, double vec[], int opt){
	// opt = 1 matriz -> vetor
	// opt = 2 vetor -> matriz
	int i;
	if (opt == 1){
		for (i = 0;i<nrows;i++){
			vec[i] = A[i][ncols];
		}
	}
	if (opt == 2){
		for (i = 0;i<nrows;i++){
			A[i][ncols] = vec[i];
		}
	}
}

// ATRIB COLUMN N
void atr_col_N_Ts(double A[][1001], int ncols, int nrows, double vec[], int opt){
	// opt = 1 matriz -> vetor
	// opt = 2 vetor -> matriz
	int i;
	if (opt == 1){
		for (i = 0;i<nrows;i++){
			vec[i] = A[i][ncols];
		}
	}
	if (opt == 2){
		for (i = 0;i<nrows;i++){
			A[i][ncols] = vec[i];
		}
	}
}

// ATRIB COLUMN 2n+1
void atr_col_2np1(double A[][15], int ncols, int nrows, double vec[], int opt){
	// opt = 1 matriz -> vetor
	// opt = 2 vetor -> matriz
	int i;
	if (opt == 1){
		for (i = 0;i<nrows;i++){
			vec[i] = A[i][ncols];
		}
	}
	if (opt == 2){
		for (i = 0;i<nrows;i++){
			A[i][ncols] = vec[i];
		}
	}
}

// MULTIPLICAÇÃO MATRIZ-MATRIZ
void mm_mult(double A[][7], double B[][7], int colsA, int rowsA, int colsB, double C[][7]){
	double Ck[rowsA];
	double Bk[colsA];
	int i;
	for (i = 0;i<colsB;i++){
		atr_col(B,i,colsA,Bk,1);
		mv_mult(A,Bk,rowsA,colsA,Ck);
		atr_col(C,i,rowsA,Ck,2);
	}
}

// MULTIPLICAÇÃO MATRIZ-MATRIZ me x n
void mm_mult_me_n(double A[][4], double B[][7], int colsA, int rowsA, int colsB, double C[][7]){
	double Ck[rowsA];
	double Bk[colsA];
	int i;
	for (i = 0;i<colsB;i++){
		atr_col(B,i,colsA,Bk,1);
		mv_mult_me(A,Bk,rowsA,colsA,Ck);
		atr_col(C,i,rowsA,Ck,2);
	}
}

// MULTIPLICAÇÃO MATRIZ-MATRIZ ME
void mm_mult_me(double A[][4], double B[][4], int colsA, int rowsA, int colsB, double C[][4]){
	double Ck[rowsA];
	double Bk[colsA];
	int i;
	for (i = 0;i<colsB;i++){
		atr_col_me(B,i,colsA,Bk,1);
		mv_mult_me(A,Bk,rowsA,colsA,Ck);
		atr_col_me(C,i,rowsA,Ck,2);
	}
}

// MULTIPLICAÇÃO MATRIZ-MATRIZ 2n+1
void mm_mult_2np1(double A[][15], double B[][15], int colsA, int rowsA, int colsB, double C[][15]){
	double Ck[rowsA];
	double Bk[colsA];
	int i;
	for (i = 0;i<colsB;i++){
		atr_col_2np1(B,i,colsA,Bk,1);
		mv_mult_2np1(A,Bk,rowsA,colsA,Ck);
		atr_col_2np1(C,i,rowsA,Ck,2);
	}
}

// MULTIPLICAÇÃO MATRIZ-MATRIZ 2n+1b
void mm_mult_2np1b(double A[][15], double B[][7], int colsA, int rowsA, int colsB, double C[][7]){
	double Ck[rowsA];
	double Bk[colsA];
	int i;
	for (i = 0;i<colsB;i++){
		atr_col(B,i,colsA,Bk,1);
		mv_mult_2np1(A,Bk,rowsA,colsA,Ck);
		atr_col(C,i,rowsA,Ck,2);
	}
}

// MULTIPLICAÇÃO MATRIZ-MATRIZ 2n+1b
void mm_mult_2np1b_me(double A[][15], double B[][4], int colsA, int rowsA, int colsB, double C[][4]){
	double Ck[rowsA];
	double Bk[colsA];
	int i;
	for (i = 0;i<colsB;i++){
		atr_col_me(B,i,colsA,Bk,1);
		mv_mult_2np1(A,Bk,rowsA,colsA,Ck);
		atr_col_me(C,i,rowsA,Ck,2);
	}
}

// DETERMINANTE
double determinant(double a[][7],int k)
{
  double s=1,det=0.0,b[n][n];
  int i,j,m,k1,c;
  if (k==1)
    {
     return (a[0][0]);
    }
  else
    {
     det=0;
     for (c=0;c<k;c++)
       {
        m=0;
        k1=0;
        for (i=0;i<k;i++)
          {
            for (j=0;j<k;j++)
              {
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][k1]=a[i][j];
                   if (k1<(k-2))
                    k1++;
                   else
                    {
                     k1=0;
                     m++;
                     }
                   }
               }
             }
          det=det + s * (a[0][c] * determinant(b,k-1));
          s=-1 * s;
          }
    }
 
    return (det);
}

// DETERMINANTE ME
double determinant_me(double a[][4],int k)
{
  double s=1,det=0.0,b[n_me][n_me];
  int i,j,m,k1,c;
  if (k==1)
    {
     return (a[0][0]);
    }
  else
    {
     det=0;
     for (c=0;c<k;c++)
       {
        m=0;
        k1=0;
        for (i=0;i<k;i++)
          {
            for (j=0;j<k;j++)
              {
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][k1]=a[i][j];
                   if (k1<(k-2))
                    k1++;
                   else
                    {
                     k1=0;
                     m++;
                     }
                   }
               }
             }
          det=det + s * (a[0][c] * determinant_me(b,k-1));
          s=-1 * s;
          }
    }
 
    return (det);
}

// COFATOR
void cofactor(double num[][7],double f,double fac[][7])
{
 double b[n][n];
 int p,q,m,k1,i,j;
 for (q=0;q<f;q++)
 {
   for (p=0;p<f;p++)
    {
     m=0;
     k1=0;
     for (i=0;i<f;i++)
     {
       for (j=0;j<f;j++)
        {
          if (i != q && j != p)
          {
            b[m][k1]=num[i][j];
            if (k1<(f-2))
             k1++;
            else
             {
               k1=0;
               m++;
               }
            }
        }
      }
      fac[q][p]=pow(-1,q + p) * determinant(b,f-1);
    }
  }
  
}

// COFATOR ME
void cofactor_me(double num[][4],double f,double fac[][4])
{
 double b[n_me][n_me];
 int p,q,m,k1,i,j;
 for (q=0;q<f;q++)
 {
   for (p=0;p<f;p++)
    {
     m=0;
     k1=0;
     for (i=0;i<f;i++)
     {
       for (j=0;j<f;j++)
        {
          if (i != q && j != p)
          {
            b[m][k1]=num[i][j];
            if (k1<(f-2))
             k1++;
            else
             {
               k1=0;
               m++;
               }
            }
        }
      }
      fac[q][p]=pow(-1,q + p) * determinant_me(b,f-1);
    }
  }
  
}

// MATRIZ INVERSA
void inverse(double A[][7], int r, double invA[][7])
{
  int i,j;
  double b[n][n],d;
  double fac[n][n];
  cofactor(A,n,fac); 
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
         b[i][j]=fac[j][i];
        }
    }
  d=determinant(A,r);
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
        invA[i][j]=b[i][j] / d;
        }
    }

}

// MATRIZ INVERSA ME
void inverse_me(double A[][4], int r, double invA[][4])
{
  int i,j;
  double b[n_me][n_me],d;
  double fac[n_me][n_me];
  cofactor_me(A,n_me,fac); 
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
         b[i][j]=fac[j][i];
        }
    }
  d=determinant_me(A,r);
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
        invA[i][j]=b[i][j] / d;
        }
    }

}

// MEASURING FUNCTION
void measf (double x[], double h1[])
{   
      h1[0] = x[0];
      h1[1] = x[1];
      h1[2] = x[2];
      h1[3] = x[5];
}

// TRANSPOSE
void transp(double A[][7], double At[][7], int nrows,int ncols){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			At[i][j] = A[j][i];
		}
	}
}

// TRANSPOSE n_me
void transp_2np1_me (double A[][15], double At[][4], int nrows,int ncols){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			At[i][j] = A[j][i];
		}
	}
}

// TRANSPOSE 2np1
void transp_2np1(double A[][15], double At[][7], int nrows,int ncols){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			At[i][j] = A[j][i];
		}
	}
}

// TRANSPOSE 2np1
void transp_me(double A[][4], double At[][4], int nrows,int ncols){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			At[i][j] = A[j][i];
		}
	}
}

// TRANSPOSE n_me
void transp_n_me(double A[][4], double At[][7], int nrows,int ncols){
	int i,j;
	for (i = 0;i<nrows;i++){
		for (j = 0;j<ncols;j++){
			At[i][j] = A[j][i];
		}
	}
}

// CHOLESKY DECOMPOSITION
void cholesky (double A[][7], double A_chol[][7], int n1)
{
   int i, j, k;
   int retval = 1;
   double At[n][n];
   double buffA[n][n];
  
   transp(A,At,n,n);
   mm_add(A,At,n,n,buffA);
   me_mult(buffA,n,n,0.5);
   
   for (i = 0; i < n1; i++) {
      A_chol[i][i] = buffA[i][i];
      for (k = 0; k<i; k++){
	 	A_chol[i][i] = A_chol[i][i] - A_chol[k][i]*A_chol[k][i];
      	}
	
	A_chol[i][i] = sqrt(A_chol[i][i]);
	
     for (j = i+1; j<n1; j++) {
	 	A_chol[i][j] = A[i][j];
	 		for (k = 0; k<i; k++){
	    		A_chol[i][j] = A_chol[i][j] - A_chol[k][i]*A_chol[k][j];
      		}		
      	A_chol[i][j] = A_chol[i][j]/A_chol[i][i];
      	
	}
	//cout << endl;
	}
	//system("pause");
	
	for (i=1;i<n1;i++){
		for (j=0;j<i;j++){
			A_chol[i][j] = 0.0;
		}
	}
}

// ESTIMATED STATES
void statef_est(double x[],double Ts,double u[],double f[]){
	// Definição dos parâmetros
	double Rs = 0.1;
	//double Rs = x[7];
	double Rfd = 0.016;
	//double Rfd = x[8];
	double Rkd = 0.17;
	double Rkq = 0.17;
	double Lls = 0.00079;
	//double Lmd = 0.002;
	double Lmq = 0.002;
	double Lmd = 0.002;
	//double Lmq = x[7];
	double Llfd = 0.00037;
	double Llkd = 0.00028;
	double Llkq = 0.00091;
	double J = 0.2;
	double B = 0.01;
	//double B = x[8];
	double buff[n];
	double buff2[n];
	double buff3[n];
	
	double L_dq_inv[7][7] = {{- 1/(Lls + Lmq) - pow(Lmq,2)/(pow(Lls + Lmq,2)*(Llkq + Lmq - pow(Lmq,2)/(Lls + Lmq))),0, Lmq/((Lls + Lmq)*(Llkq + Lmq - pow(Lmq,2)/(Lls + Lmq))),0,0, 0, 0},
               {0,- pow((Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))),2)/(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))) - 1/(Lls + Lmd) - pow(Lmd,2)/(pow((Lls + Lmd),2)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, Lmd/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))) - ((Lmd - pow(Lmd,2)/(Lls + Lmd))*(Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))))/((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))), (Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))))/(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, 0},
               {-Lmq/((Lls + Lmq)*(Llkq + Lmq - pow(Lmq,2)/(Lls + Lmq))),0,1/(Llkq + Lmq - pow(Lmq,2)/(Lls + Lmq)),0,0, 0, 0},
            {0, ((Lmd - pow(Lmd,2)/(Lls + Lmd))*(Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))))/((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))) - Lmd/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, 1/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)) + pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(pow((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)),2)*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))), -(Lmd - pow(Lmd,2)/(Lls + Lmd))/((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))), 0, 0},
            {0,-(Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))))/(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, -(Lmd - pow(Lmd,2)/(Lls + Lmd))/((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))), 1/(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, 0},
            {0,0,0,0,0,1.,0},
            {0,0,0,0,0,0,1.}};
            

            
    double R_dq[7][7] = {{Rs, 0.5*Poles*x[5]*(Lls+Lmd),0, -0.5*Poles*x[5]*Lmd, -0.5*Poles*x[5]*Lmd,0,0},
            {-0.5*Poles*x[5]*(Lls+Lmq),Rs,0.5*Poles*x[5]*Lmq, 0, 0, 0, 0},
         	{0, 0, -Rkq, 0, 0, 0, 0},
            {0, 0, 0, -Rkd, 0, 0, 0},
            {0, 0, 0, 0, -Rfd, 0, 0},
            {0, 0, 0, 0, 0, -B/J, 0},
            {0, 0, 0, 0, 0, 0.5*Poles,0}};
            
//        	cout << "R_dq" << endl << endl;
//            for (int i=0;i<n;i++){
//            	for (int j=0;j<n;j++){
//            		cout << R_dq[i][j] << "  " ;
//				}
//				cout << endl;
//			}
//			cout << endl << endl;
//			system("pause");
//			system("cls");
            
	double Te = (3/2)*(Poles/2)*(Lmd*(-x[1]+x[3]+x[4])*x[0] - Lmq*(-x[0]+x[2])*x[1]);  
	mv_mult(R_dq,x,n,n,buff);
	add(buff,u,n,buff);
	mv_mult(L_dq_inv,buff,n,n,buff2);
	sc_vec(Ts,buff2,n,buff3);
	add(x,buff3,n,f);
	//for (int i=0;i<n;i++){
	//	cout << f[i] << "  ";
	//}
	//cout << endl << endl;
	
}

// STATES
//double statef(double x[],double Ts,double u[],double f[]){
	// Definição dos parâmetros
//	double Rs = 0.1;
//	double Rfd = 0.016;
//	double Rkd = 0.17;
//	double Rkq = 0.17;
//	double Lls = 0.00079;
//	double Lmd = 0.002;
//	double Lmq = 0.002;
//	double Llfd = 0.00037;
//	double Llkq = 0.00028;
//	double Llkd = 0.00091;
//	double J = 0.2;
//	double B = 0.01;
	//double L_dq_inv[n_orig][n_orig];
//	double buff[n_orig];
//	double buff2[n_orig];
//	double buff3[n_orig];
	
//	double L_dq[n_orig][n_orig] = {{-(Lls+Lmq), 0, Lmq, 0, 0, 0, 0},
//       {0, -(Lls+Lmd), 0, Lmd, Lmd, 0, 0},
//       {-Lmq, 0, (Llkq+Lmq), 0, 0, 0, 0},
//       {0, -Lmd, 0, (Llkd+Lmd), Lmd, 0, 0},
//       {0, -Lmd, 0, Lmd, (Llfd+Lmd), 0, 0},
//       {0, 0, 0, 0, 0, 1., 0},
//       {0, 0, 0, 0, 0, 0, 1.}};

//	double L_dq_inv[7][7] = {{- 1/(Lls + Lmq) - pow(Lmq,2)/(pow(Lls + Lmq,2)*(Llkq + Lmq - pow(Lmq,2)/(Lls + Lmq))),0, Lmq/((Lls + Lmq)*(Llkq + Lmq - pow(Lmq,2)/(Lls + Lmq))),0,0, 0, 0},
//              {0,- pow((Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))),2)/(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))) - 1/(Lls + Lmd) - pow(Lmd,2)/(pow((Lls + Lmd),2)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, Lmd/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))) - ((Lmd - pow(Lmd,2)/(Lls + Lmd))*(Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))))/((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))), (Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))))/(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, 0},
//               {-Lmq/((Lls + Lmq)*(Llkq + Lmq - pow(Lmq,2)/(Lls + Lmq))),0,1/(Llkq + Lmq - pow(Lmq,2)/(Lls + Lmq)),0,0, 0, 0},
//            {0, ((Lmd - pow(Lmd,2)/(Lls + Lmd))*(Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))))/((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))) - Lmd/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, 1/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)) + pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(pow((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)),2)*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))), -(Lmd - pow(Lmd,2)/(Lls + Lmd))/((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))), 0, 0},
////            {0,-(Lmd/(Lls + Lmd) - (Lmd*(Lmd - pow(Lmd,2)/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))))/(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, -(Lmd - pow(Lmd,2)/(Lls + Lmd))/((Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))*(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd)))), 1/(Llfd + Lmd - pow(Lmd,2)/(Lls + Lmd) - pow((Lmd - pow(Lmd,2)/(Lls + Lmd)),2)/(Llkd + Lmd - pow(Lmd,2)/(Lls + Lmd))), 0, 0},
//            {0,0,0,0,0,1.,0},
//            {0,0,0,0,0,0,1.}};
       
///	double R_dq[7][7] = {{Rs, 0.5*Poles*x[5]*(Lls+Lmd), 0, -0.5*Poles*x[5]*Lmd, -0.5*Poles*x[5]*Lmd, 0, 0},
//       {-0.5*Poles*x[5]*(Lls+Lmq),Rs, 0.5*Poles*x[5]*Lmq, 0, 0, 0, 0},
//        {0, 0, -Rkq, 0, 0, 0, 0},
//        {0, 0, 0, -Rkd, 0, 0, 0},
//        {0, 0, 0, 0, -Rfd, 0, 0},
//        {0, 0, 0, 0, 0, -B/J, 0},
//        {0, 0, 0, 0, 0, 0.5*Poles, 0}};
       
//	double Te = (3/2)*(Poles/2)*(Lmd*(-x[1]+x[3]+x[4])*x[0] - Lmq*(-x[0]+x[2])*x[1]);  

//	mv_mult_orig(R_dq,x,n_orig,n_orig,buff);
//	add(buff,u,n_orig,buff);
//	//inverse(L_dq,n_orig,L_dq_inv);
//	mv_mult_orig(L_dq_inv,buff,n_orig,n_orig,buff2);
//	sc_vec(Ts,buff2,n_orig,buff3);
//	add(x,buff3,n_orig,f);
	
	//  for (int i=0;i<n;i++){
    //  	for (int j=0;j<n;j++){
    //  		cout << L_dq_inv[i][j] << "    " ;
	//	   }
	//	    cout << endl;
	//  }
//}
void me_mult3(double A[][3],int rows, int cols, double c){
	int i,j;
	for (i = 0;i<rows;i++){
		for (j = 0;j<cols;j++){
			A[i][j] = c*A[i][j];
		}
	}
}

void mv_mult3(double A[][3], double x[], int rows, int cols, double b[]){ 
	int i,j;
    for (i = 0; i < rows; i++)
    {
    	b[i] = 0.0;
    	for (j = 0;j<cols;j++){
			b[i] += (A[i][j] * x[j]);
			}
			
    } 
}

void abc2dqn(double Vabc[],double theta,double Vdqn[]){
	//double pi = 3.141592;
	double Mdqn[3][3] = {{sin(theta),sin(theta-(2.0/3.0)*pi),sin(theta+(2.0/3.0)*pi)},{cos(theta),cos(theta-(2.0/3.0)*pi),cos(theta+(2.0/3.0)*pi)},{sqrt(2)/2, sqrt(2)/2, sqrt(2)/2}};
	me_mult3(Mdqn,3,3,(2.0/3.0));
	//cout << Mdqn[0][0] << " "	<< Mdqn[0][1] << " " << Mdqn[0][2] << endl;
	//cout << Mdqn[1][0] << " "	<< Mdqn[1][1] << " " << Mdqn[1][2] << endl;
	//cout << Mdqn[2][0] << " "	<< Mdqn[2][1] << " " << Mdqn[2][2] << endl;
	mv_mult3(Mdqn,Vabc,3,3,Vdqn);
}

void dqn2abc(double Vdqn[],double theta,double Vabc[]){
	//double pi = 3.141592;
	double Mabc[3][3] = {{sin(theta),cos(theta), sqrt(2)/2},{sin(theta-(2.0/3.0)*pi),cos(theta-(2.0/3.0)*pi),sqrt(2)/2},{sin(theta+(2.0/3.0)*pi),cos(theta+(2.0/3.0)*pi),sqrt(2)/2}};
	//cout << Mabc[0][0] << " "	<< Mabc[0][1] << " " << Mabc[0][2] << endl;
	//cout << Mabc[1][0] << " "	<< Mabc[1][1] << " " << Mabc[1][2] << endl;
	//cout << Mabc[2][0] << " "	<< Mabc[2][1] << " " << Mabc[2][2] << endl;
	mv_mult3(Mabc,Vdqn,3,3,Vabc);
}

#endif // HEADER_H_

