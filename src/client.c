/**** 
	FaRSA

	version: 2.0

	Authors:

	Tianyi Chen, Frank E. Curtis, Daniel P. Robinson
	
	April, 2018
*/

#include "farsa.h"

int main( int argc, char **argv )
{
	/* input structure of FaRSA */
	struct Input_FaRSA input_farsa;

	/* input structure of FaRSA */
	struct Output_FaRSA output_farsa;

	/* run farsa */
	output_farsa 	= farsa( argc, argv, &input_farsa );
	
	/* 0 is optimal solution is found, 1 is reach max_iter, 2 is reach max_time */
	return output_farsa.status;
	
}