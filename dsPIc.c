#include <stdlib.h>
#include<math.h>
//#include<iostream>
#include<stdio.h>

//using namespace std;

static double tempo = 0;
static double t = 0;
static float Va;
static float Vb;
static float Ia;
static float Ib;
static float Vf;
static float theta_e;

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
//float in[7] = {1,1,1,1,1,1,1};
//float out[12];
//float delt = 0.00001;
int flag = 0;

// DEC2BIN
//void dec2bin(int decimal,int A[]){
//	int i, k;
//	
//// posição 0
//for (i = n_bits-1;i>=0;i--){
//
//	k = decimal >> i;
//	if (k & 1){
//		A[i] = 1;
//	}
//	else{
//		A[i] = 0;
//	}
//}
//}
//
//// BIN2DEC
//int bin2dec(int A[]){
//	int i;
//	int decimal = 0;
//// posição 0
//for (i = 0;i<n_bits;i++){
//	decimal = decimal + pow(2,i)*A[i];
//	}
//	return decimal;
//}

int main(){
	while (1){

tempo = tempo + delt;
t = t + delt;

if (tempo>0.0001){
	flag = 1;
	
//contador
int i,j;

// Leitura dos dados
Va = round(in[0]/ampl_Va*(pow(2,n_bits)-1));
Vb = round(in[1]/ampl_Vb*(pow(2,n_bits)-1));
Ia = round(in[2]/ampl_Ia*(pow(2,n_bits)-1));
Ib = round(in[3]/ampl_Ib*(pow(2,n_bits)-1));
Vf = round(in[4]/ampl_Vf*(pow(2,n_bits)-1));
theta_e = round(in[5]/ampl_theta_e*(pow(2,n_bits)-1));

// INICIO //
while (in[6]==0){
	C[0] = 1;
	C[1] = 1;
	}
///////////

// VA ///////////////////////////////////////////////////////// VA ///////////////////////////////////////////
while (in[6]==1){

	dec2bin(Va,A);
	// seleção do parametro
	C[0] = 1;
	C[1] = 0;


	for (i=0;i<n_bits;i++){
		out[i] = A[i];
	}

	for (i=n_bits;i<n_bits+2;i++){
		out[i] = C[i-n_bits];
	}
	//for (i=0;i<n_bits+3;i++){
	//	cout << out[i] << endl;
	//}
}
/// FIM VA /////////////////////////////////////////////////////// FIM VA ////////////////////////

// VB///////////////////////////////////////////////////////// VB ///////////////////////////////////////////
while (in[6]==0){

	dec2bin(Vb,A);
	// seleção do parametro
	C[0] = 0;
	C[1] = 1;

	for (i=0;i<n_bits;i++){
		out[i] = A[i];
	}

	for (i=n_bits;i<n_bits+2;i++){
		out[i] = C[i-n_bits];
	}

	//for (i=0;i<n_bits+3;i++){
	//	cout << out[i] << endl;
	//}
}
/// FIM VB /////////////////////////////////////////////////////// FIM VB ////////////////////////

// IA ///////////////////////////////////////////////////////// IA ///////////////////////////////////////////
while (in[6]==1){

	dec2bin(Ia,A);
	// seleção do parametro
	C[0] = 1;
	C[1] = 0;

	for (i=0;i<n_bits;i++){
		out[i] = A[i];
	}

	for (i=n_bits;i<n_bits+2;i++){
		out[i] = C[i-n_bits];
	}
	//for (i=0;i<n_bits+3;i++){
	//	cout << out[i] << endl;
	//}
}
/// FIM IA /////////////////////////////////////////////////////// FIM IA ////////////////////////

// IB ///////////////////////////////////////////////////////// IB ///////////////////////////////////////////
while (in[6]==0){

	dec2bin(Ib,A);
	// seleção do parametro
	C[0] = 0;
	C[1] = 1;

	for (i=0;i<n_bits;i++){
		out[i] = A[i];
	}

	for (i=n_bits;i<n_bits+2;i++){
		out[i] = C[i-n_bits];
	}
	//for (i=0;i<n_bits+3;i++){
	//	cout << out[i] << endl;
	//}
}
/// FIM IB /////////////////////////////////////////////////////// FIM IB ////////////////////////

// Vf ///////////////////////////////////////////////////////// Vf ///////////////////////////////////////////
while (in[6]==1){

	dec2bin(Vf,A);
	// seleção do parametro
	C[0] = 1;
	C[1] = 0;

	for (i=0;i<n_bits;i++){
		out[i] = A[i];
	}

	for (i=n_bits;i<n_bits+2;i++){
		out[i] = C[i-n_bits];	
	}
	//for (i=0;i<n_bits+3;i++){
	//	cout << out[i] << endl;
	//}
}
/// FIM Vf /////////////////////////////////////////////////////// FIM Vf ////////////////////////

// theta ///////////////////////////////////////////////////////// theta ///////////////////////////////////////////
while (in[6]==0){

	dec2bin(theta_e,A);
	// seleção do parametro
	C[0] = 0;
	C[1] = 1;

	for (i=0;i<n_bits;i++){
		out[i] = A[i];
	}

	for (i=n_bits;i<n_bits+2;i++){
		out[i] = C[i-n_bits];
	}
	//for (i=0;i<n_bits+3;i++){
	//	cout << out[i] << endl;
	//}
}
/// FIM theta /////////////////////////////////////////////////////// FIM theta ////////////////////////
while (in[6]==1){

	C[0] = 0;
	C[1] = 0;
	flag = 0;
}
}
}
}
