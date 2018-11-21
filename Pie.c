#include <stdlib.h>
#include<math.h>
//#include<iostream>
#include<stdio.h>

using namespace std;

static float Va;
static float Vb;
static float Vc;
static float Ia;
static float Ib;
static float Ic;
static float Vf;
static float If;
static float theta_e;
static float Vabc[3];
static float Vdq0[3];
static float Iabc[3];
static float Idq0[3];
const float Rfd = 0.016;

// binary
static int A[10];
static int C[2];

// dados
const int ampl_Va = 311;
const int ampl_Vb = 311;
const int ampl_Ia = 1;
const int ampl_Ib = 1;
const int ampl_Vf = 10;
const int ampl_theta_e = 2*3.1415926;
const int n_bits = 10;
//float in[12] = {1,1,1,1,1,1,1};
//float out[8];
int flag = 0;

// DEC2BIN
void dec2bin(int decimal,int A[]){
	int i, k;
	
// posição 0
for (i = n_bits-1;i>=0;i--){

	k = decimal >> i;
	if (k & 1){
		A[i] = 1;
	}
	else{
		A[i] = 0;
	}
}
}

// BIN2DEC
int bin2dec(int A[]){
	int i;
	int decimal = 0;
// posição 0
for (i = 0;i<n_bits;i++){
	decimal = decimal + pow(2,i)*A[i];
	}
	return decimal;
}

int main(){
	while (1){

//contador
int i,j;

// Leitura dos dados
//Va = round(in[0]/ampl_Va*(pow(2,n_bits)-1));
//Vb = round(in[1]/ampl_Vb*(pow(2,n_bits)-1));
//Ia = round(in[2]/ampl_Ia*(pow(2,n_bits)-1));
//Ib = round(in[3]/ampl_Ib*(pow(2,n_bits)-1));
//Vf = round(in[4]/ampl_Vf*(pow(2,n_bits)-1));
//theta_e = round(in[5]/ampl_theta_e*(pow(2,n_bits)-1));

// INICIO //

///////////
while (out[7]==0){
	C[0] = in[10];
	C[1] = in[11];
	if (C[0]==1 && C[1]==1){
		out[7] = 1;
	}
}
// VA ///////////////////////////////////////////////////////// VA ///////////////////////////////////////////
while (out[7]==1){
	C[0] = in[10];
	C[1] = in[11];
	if (C[0]==1 && C[1]==0){
	for (i=0;i<n_bits;i++){
		A[i] = in[i];
	}
	}
	Va = bin2dec(A)/(pow(2,n_bits)-1)*ampl_Va;
	out[7] = 0;
}
/// FIM VA /////////////////////////////////////////////////////// FIM VA ////////////////////////

// VB///////////////////////////////////////////////////////// VB ///////////////////////////////////////////
while (out[7]==0){
	C[0] = in[10];
	C[1] = in[11];
	if (C[0]==0 && C[1]==1){
	for (i=0;i<n_bits;i++){
		A[i] = in[i];
	}
	}
	Vb = bin2dec(A)/(pow(2,n_bits)-1)*ampl_Vb;
	out[7] = 1;

}
/// FIM VB /////////////////////////////////////////////////////// FIM VB ////////////////////////

// IA ///////////////////////////////////////////////////////// IA ///////////////////////////////////////////
while (out[7]==1){
	C[0] = in[10];
	C[1] = in[11];
	if (C[0]==1 && C[1]==0){
	for (i=0;i<n_bits;i++){
		A[i] = in[i];
	}
	}
	Ia = bin2dec(A)/(pow(2,n_bits)-1)*ampl_Ia;
	out[7] = 0;

}
/// FIM IA /////////////////////////////////////////////////////// FIM IA ////////////////////////

// IB ///////////////////////////////////////////////////////// IB ///////////////////////////////////////////
while (out[7]==0){
	C[0] = in[10];
	C[1] = in[11];
	if (C[0]==0 && C[1]==1){
	for (i=0;i<n_bits;i++){
		A[i] = in[i];
	}
	}
	Ib = bin2dec(A)/(pow(2,n_bits)-1)*ampl_Ib;
	out[7] = 1;

}
/// FIM IB /////////////////////////////////////////////////////// FIM IB ////////////////////////

// Vf ///////////////////////////////////////////////////////// Vf ///////////////////////////////////////////
while (out[7]==1){
	C[0] = in[10];
	C[1] = in[11];
	if (C[0]==1 && C[1]==0){
	for (i=0;i<n_bits;i++){
		A[i] = in[i];
	}
	}
	Vf = bin2dec(A)/(pow(2,n_bits)-1)*ampl_Vf;
	out[7] = 0;

}
/// FIM Vf /////////////////////////////////////////////////////// FIM Vf ////////////////////////

// theta ///////////////////////////////////////////////////////// theta ///////////////////////////////////////////
while (out[7]==0){
	C[0] = in[10];
	C[1] = in[11];
	if (C[0]==0 && C[1]==1){
	for (i=0;i<n_bits;i++){
		A[i] = in[i];
	}
	}
	theta_e = bin2dec(A)/(pow(2,n_bits)-1)*ampl_theta_e;
	out[7] = 1;

}
/// FIM theta /////////////////////////////////////////////////////// FIM theta ////////////////////////
while (out[7]==1){
	C[0] = in[10];
	C[1] = in[11];
	if (C[0]==0 && C[1]==1){
	out[7] = 0;
	flag = 1;
}
}
if (flag == 1){
	//for (k=0;k<N;k++){
		
//		s[5] = -2.0*3.141592*60/(4.0*0.5);

		Vc = -Va - Vb;
		Ic = -Ia - Ic;
		// montagem Vabc e Iabc
		Vabc[0] = Va;
		Vabc[1] = Vb;
		Vabc[2] = Vc;
		Iabc[0] = Ia;
		Iabc[1] = Ib;
		Iabc[2] = Ic;
		If = Vf/Rfd;
		
		abc2dqn(Vabc,theta_e,Vdq0);
		abc2dqn(Iabc,theta_e,Idq0);
		
		// entrada
		u[0] = Vdq0[1];
		u[1] = Vdq0[0];
		u[4] = Vf;
		
		// estados
		s[0] = Idq0[1];
		s[1] = Idq0[0];
		s[2] = -2*pi*f;
		rdn_vec(qn,n_orig);
		sc_vec(q,qn,n_orig,buff);
		add(s,buff,n_orig,s);
		zeros_vec(buff,n_orig);
		
		measf(s,z);	
		rdn_vec(rn,n_me);
		sc_vec(r,rn,n_me,buff1);
		add(z,buff1,n_me,z);

		zeros_vec(buff1,n_me);
		

// UFK ****************************************************

	Wm[0] = lambda/c;
	double c1;
	Wc[0] = Wm[0] + (1 - pow(alpha,2)+beta);
	for (i = 1; i<2*n+1; i++){
		Wm[i] = 0.5/c;
		Wc[i] = 0.5/c;
	}
	c1 = sqrt(c);
	
	cholesky(P, P_chol, n);
	me_mult(P_chol,n,n,c1);
	transp(P_chol, At, n,n);

	// coluna igual a x
	for (i = 0;i<n;i++){
		X[i][0] = x[i];
	}

	// outras colunas somadas a A
	for (j = 1;j<n+1;j++){
		for (i = 0;i<n;i++){
			X[i][j] = x[i] + At[i][j-1];
		}
	}
	
	for (j = n+1;j<2*n+1;j++){
		for (i = 0;i<n;i++){
			X[i][j] = x[i] - At[i][j-(n+1)];
		}
	}
	
	// Show X
	
	zeros_vec(x1,n);
	int k1;
	for (k1 = 0;k1 <2*n+1;k1++){
			atr_col_2np1(X,k1,n,Xxk,1);
			statef_est(Xxk,Ts,u,f);
			atr_col_2np1(X1,k1,n,f,2);
			sc_vec(Wm[k1],f,n,buff8);
			add(x1,buff8,n,x1);
			zeros_vec(buff8,n);
		}

	// Y1 = X1 - x1(:,ones(1,L);
	for (k1 = 0;k1<2*n+1;k1++){
		atr_col_2np1(X1,k1,n,Yxk,1);
		subtr(Yxk,x1,n,Yxk);
		atr_col_2np1(X2,k1,n,Yxk,2);
	}

	// P = Y1*diag(Wc)*Y1' + Rn;
	
	for (i = 0;i<2*n+1;i++){
		for (j = 0;j<2*n+1;j++){
			if (i==j){
				diagWc[i][j] = Wc[i];
			}
			if (i!=j){
				diagWc[i][j] = 0.0;
			}
		}
	}
	
	mm_mult_2np1(X2,diagWc,2*n+1,n,2*n+1,Pbuff);
	transp_2np1(X2,X2t,2*n+1,n);
	mm_mult_2np1b(Pbuff,X2t,2*n+1,n,n,Px);
	mm_add(Px,Qn,n,n,Px);

// ut_meas	
	zeros_vec(z1,n_me);
	for (k1 = 0;k1 <2*n+1;k1++){
			atr_col_2np1(X1,k1,n,Xzk,1);
			measf(Xzk,h1);
			atr_col_2np1(Z1,k1,n_me,h1,2);
			sc_vec(Wm[k1],h1,n_me,buff9);
			add(z1,buff9,n_me,z1);
			zeros_vec(buff9,n_me);
		}
		
	// Y' = Y - z1(:,ones(1,L);
	for (k1 = 0;k1<2*n+1;k1++){
		atr_col_2np1(Z1,k1,n_me,Yzk,1);
		subtr(Yzk,z1,n_me,Yzk);
		atr_col_2np1(Z2,k1,n_me,Yzk,2);
	}
	
	
	// P = Z2*diag(Wc)*Z2' + Rn;
	
	for (i = 0;i<2*n+1;i++){
		for (j = 0;j<2*n+1;j++){
			if (i==j){
				diagWc[i][j] = Wc[i];
			}
			if (i!=j){
				diagWc[i][j] = 0.0;
			}
		}
	}
	
	mm_mult_2np1(Z2,diagWc,2*n+1,n_me,2*n+1,Pybuff);
	transp_2np1_me(Z2,Z2t,2*n+1,n_me);
	mm_mult_2np1b_me(Pybuff,Z2t,2*n+1,n_me,n_me,Py);
	mm_add_me(Py,Rn,n_me,n_me,Py);
	
	
	for (i = 0;i<2*n+1;i++){
		for (j = 0;j<2*n+1;j++){
			if (i==j){
				diagWc[i][j] = Wc[i];
			}
			if (i!=j){
				diagWc[i][j] = 0.0;
			}
		}
	}
	
	mm_mult_2np1(X2,diagWc,2*n+1,2*n+1,2*n+1,Pxybuff);	
	transp_2np1_me(Z2,Z2t,2*n+1,n_me);
	mm_mult_2np1b_me(Pxybuff,Z2t,2*n+1,n,n_me,Pxy);
	
	// Show Pxy
	inverse_me(Py,n_me,invPy);
	
		// Show invPy
	
	mm_mult_me(Pxy,invPy,n_me,n,n_me,K);	
	subtr(z,z1,n_me,buff11);
	mv_mult_me(K,buff11,n,n_me,buff2);

	zeros_vec(buff11,n_me);
	// end ut_meas 
	
	
	// State update
	for (i=0;i<n;i++){
		x[i] = x1[i]+buff2[i];
	}

	zeros_vec(buff2,n);
	
	transp_n_me(K,Kt,n_me,n);
	mm_mult_me_n(Py,Kt,n_me,n_me,n,buff3);
	mm_mult_me_n(K,buff3,n_me,n,n,buff4);
	mm_subtr(Px,buff4,n,n,P);
	for (i=0;i<n;i++){
		Pvec[i] = P[i][i];
	}
		
	
	zeros_matrix(buff3,n_me);
	zeros_matrix(buff4,n);
	

// END UFK ****************************************************

	zeros_vec(Pvec,n);

	out[0] = x[0];
	out[1] = x[1];
	out[2] = x[2];
	out[3] = x[3];
	out[4] = x[4];
	out[5] = x[5];
	out[6] = x[6];

}
}
}

