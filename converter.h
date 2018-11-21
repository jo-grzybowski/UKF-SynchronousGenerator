#ifndef CONVERTER_H_
#define CONVERTER_H_
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int n_bits = 10;

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

#endif // HEADER_H_

