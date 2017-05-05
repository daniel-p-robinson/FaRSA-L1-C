#include "farsa.h"

/* Personalized */
double personalized_func( double *x );
double *personalized_grad_f( double *x );
double *personalized_hessVec( double *v);

struct Input_FaRSA input_farsa;

double **D;
double *DTData;
double **DTD;
double *Dx;
double *Data;

double *ones;

int main( int argc, char **argv ){

	int n;
	int i, j;
	n = 100;
	Data = Malloc( double, n * n );
	DTData = Malloc( double, n * n );
	D = Malloc( double *, n );
	DTD =Malloc( double *, n );
	ones = Malloc( double, n );

	for( i = 0; i < n; i++ ){
		D[i] = &Data[i*n];
		DTD[i] = &DTData[i*n];
		ones[i] = 1.0;
		for( j = 0; j < n; j++ ){
			Data[i*n+j] = 1.0 / (double) n;
			DTData[i*n+j] = 3.0 / (double) n;
			if( i == j ){
				Data[i*n+j] += 1.0;
				DTData[i*n+j] += 1.0;
			}	
		}
	}

	input_farsa.n = n;
	input_farsa.func = &personalized_func;
	input_farsa.grad_f = &personalized_grad_f;
	input_farsa.hessVec = &personalized_hessVec;


	// struct Output_FaRSA output_farsa = farsa( argc, argv, NULL );

	struct Output_FaRSA output_farsa = farsa( argc, argv, &input_farsa );
	
	return 1;

}



double personalized_func( double *x ){

	int i, j;

	double func = 0.0;

	Dx = Malloc( double, n );

	for( i = 0; i < n; i++ ){

		Dx[i] = dot_n5( D[i], x, n );
	}

	func += dot_n5( Dx, Dx, n );

	func -= 4.0 * dot_n5( Dx, ones, n );

	func += lambda * l1_n5( x, n );

	return func;
}


double *personalized_grad_f( double *x ){

	int i, j;

	double *grad_f = Malloc( double, n );

	for( i = 0; i < n; i++ ){
		Dx[i] = dot_n5( D[i], x, n );
	}

	for( i = 0; i < n; i++ ){		
		grad_f[i] = 2.0 * dot_n5( D[i], Dx, n ) - 8.0;		
	}

	return grad_f;
}


double *personalized_hessVec( double *v ){
	int i, j;
	double *hv = Malloc( double, nI_phi );
	
	for( i = 0; i < nI_phi; i++ ){
		hv[i] = 2.0 * dot_n5( DTD[I_phi[i]], v, n );
	}

	return hv;
}