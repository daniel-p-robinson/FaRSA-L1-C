/* 
	FaRSA

	version: 2.0

	Authors:

	Tianyi Chen, Frank E. Curtis, Daniel P. Robinson
	
	April, 2018
*/
#ifndef farsa_h   /* Include guard */
#define farsa_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/* macro for objective function type */
#define GENERIC_FUNC 	0
#define	LOGISTIC_LOSS 	1
#define ELASTIC_NET		2
#define SMOOTH_SVM 		3

/* macro for sub-problem sover */
#define USE_DIRECT 		0
#define USE_CG 			1
#define USE_CD_QP 		2

/* macro for saving A^TA for least-square */
#define NOT_SAVE		0
#define MAY_SAVE 		1
#define MUST_SAVE 		2

#define MAX_VALUE 		1E7
#define Malloc( type, n ) ( type * ) malloc( (n) * sizeof(type) )
#define MAX( x, y ) ( ((x) > (y)) ? (x) : (y) )
#define MIN( x, y ) ( ((x) < (y)) ? (x) : (y) )
#ifndef FALSE
#define FALSE 			0
#endif
#ifndef TRUE
#define TRUE 			1
#endif

/* Sparse data node that contains the row index, column index and the value of data point*/
typedef struct Node 
{
	int 	index;
	double 	value;

} Node_t;

/* Parameter structure */
typedef struct Parameter 
{

	int 	max_iter;	
	int 	print_level;
	int 	tryCG;
	int 	tryCD;
	int 	tryCrossOver;
	int 	max_CG_iter;
	int 	max_CD_iter;
	int 	scaling;
	int 	maxback;
	int 	nVmaxAllow;
	int 	checkResEvery;
	int 	termType;
	int 	objective_function_type;
	int 	sub_solver;
	int 	print_every; /* print status every print_every iters */
	int 	do_test; /* do testing routine or not */
	int 	print_final_x; /* output final x */

	char 	train_file[1024];
	char 	profile_file[1024];
	char 	data_format[1024];
	char 	output_file[1024];
	char 	name[1024];
	char 	test_file[1024];

	double 	max_time; 
	double 	Gamma;
	double 	eta_r;
	double 	eta;
	double 	xi;
	double 	tol_absolute;
	double 	tol_relative;
	double 	betaFrac;
	double 	phiFrac;
	double 	crossOverTol;
	double 	TRradiusPhi;
	double 	TRradiusBeta;
	double 	fracViol;
	double 	lambda;
	double 	lambdaCoeff;
	double 	ttol;
	
} Parameter_t ;

/* Input structure for farsa */
typedef struct Input_FaRSA
{

	int 	n; /* n is the number of variables */
	double 	( *func )( double * ); /* func is to return objective function value */
	double 	*( *grad_f )( double * ); /* calculate gradient */
	void 	( *hessVec )( double *, double * ); /* calculate hessian vector product */

} Input_FaRSA_t;

/* Output structure for farsa */
typedef struct Output_FaRSA
{

	double 	*x_initial;
	double 	*x_final;
	double 	F_initial;
	double 	F_final;
	double 	f_final;
	double 	x_l1_norm_final;
	double 	norm_beta_initial;
	double 	norm_beta_final;
	double 	norm_phi_initial;
	double 	norm_phi_final;
	double 	error;
	double 	ttol;
	double 	run_time;
	double 	lambda;
	double 	zero_perc;
	int 	n;
	int 	iter;
	int 	beta_iter;
	int 	phi_iter;
	int 	f_eval;
	int 	grad_eval;
	int 	hessVec_eval;
	int 	term;
	int 	status; /* 0 is optimal solution is found, 1 is reach max_iter, 2 is reach max_time */

} Output_FaRSA_t;

/* Projection result */
typedef struct Proj_output
{	
	int 	same_sign;
	int 	nV;

} Proj_output_t;

/* Objective function value result */
typedef struct Func_output
{
	
	double 	f_value;
	double	F_value;
	double 	x_l1_norm;

} Func_output_t;

/* Problem structure to save information */
typedef struct Problem 
{

	/* for training */
	int 		num_features; /* number of variables */
	int 		num_samples;
	int 		*nnz_cols; /* the number of non-zero entries for each column */
	int 		*nnz_rows; /* the number of non-zero entries for each row */
	int 		nnz; /* total number of non-zero entries */
	double 		*y;
	Node_t 		**rX; 
	Node_t 		**cX;
	Node_t 		**cXS;
	Node_t 		**rXS;

	/* for testing if needed */
	double 		*predict_y;

} Problem_t;

/* structure of CG input */
typedef struct Input_CG {

	int 	*I;
	int 	*I_p;
	int 	*I_n;
	int 	nI;
	int 	nVmaxAllow;
	int 	maxCG_iter;
	int 	n;
	int 	iter;
	double 	*x;
	double 	*grad_F;
	double 	*d_full;
	double 	eta_r;
	double 	rsType;
	double 	fracViol;
	double 	TRradius;
	double 	HvProds;
	void 	( *hessVecProd )( double *, double * );

} Input_CG_t;


/* structure of input for CD sub-problem solver */
typedef struct Input_CD {

	int 	*I;
	int 	*I_p;
	int 	*I_n;
	int 	*x_I_pos;
	int 	*x_I_neg;
	int 	nI;
	int 	nVmaxAllow;
	int 	max_CD_iter;
	int 	n;
	int 	m;
	double 	*x;
	double 	*x_I;
	double 	*grad_F;
	double 	*grad_F_I;
	double 	*help_vector; /* help vector for different objective function */
	double 	eta_r;
	double 	rsType;
	double 	fracViol;
	double 	TRradius;
	double 	HvProds;
	double 	*( *hessVecProd )( double * );
	double 	regularize;

} Input_CD_t;

/* sub-problem solver output */
typedef struct Output_sp {

	int 	nV;
	int 	nVmax;
	int 	iter;
	double 	*d_full;
	double 	res_target;
	double 	res;
	double 	dir_der;
	double 	norm_d;
	char 	sub_prob_flag[1024];

} Output_sp_t;


typedef struct Output_CD {

	int 	nV;
	int 	iter;
	double 	*d_full;
	double 	dir_der;
	char 	sub_prob_flag[1024];

} Output_CD_t;


typedef struct Test_output {
	
	int num_samples;
	int tp;
	int tn;
	int fp;
	int fn;
	double accuracy;
	double precision;
	double recall;
	double f1_score;

} Test_output_t;

/* Global variables */
int 	iter;
int 	beta_iter;
int 	phi_iter;
int 	m;
int 	n;
int 	rsType; /* Target Residual Type in phi step */
int 	iter_type; /* 1 is phi step, 0, is beta step */
int 	*I_phi;
int 	*I_phi_selected; /* selected index set for sub-problem */
int 	*I_beta;
int 	*I_beta_selected;	
int 	*I_z;
int 	*I_n;
int 	*I_p;
int 	*row_ptr_sp;
int 	n5; /* n5 use for accelerate calculation process */
int 	nI_phi;
int 	nI_phi_selected; /* selected sub-problem size */
int 	nI_beta;
int 	nI_beta_selected; 
int 	nnz; /* counter of nonzero entries */
int 	f_eval;
int 	grad_eval;
int 	hessVec_eval;
int 	use_CG;
int 	use_direct;
int 	use_CD_qp;
int 	sub_iters;

double 	*beta;
double 	*beta_I;
double 	*Beta; 
double 	*phi;
double 	*phi_I;
double 	*grad_f;
double 	*grad_F;
double 	norm_beta;
double 	norm_Beta;
double 	norm_phi;
double 	norm_beta0;
double 	norm_phi0;
double 	lambda;
double 	Gamma;
double 	*x;
double 	ttol;
double 	tol_relative;
double 	tol_absolute;
double 	alpha_B;
double 	m_d; /* double type of m */
double 	n_d; /* double type of n */

/* saving the diagonal entries of D in X^TDX for logistic loss */
double 	*diagonal;
double 	*y;
double 	*y_B;
double 	*d;
double 	*help_hv_logis; /* help vector for logistic cost*/
double 	TRradiusBeta;
double 	TRradiusPhi;
double 	alpha;
double 	*x_linesearch;
double 	dirDer;
double 	F_old;  /* objective function value for last step */
double 	F; /* objective function value on current step */
double 	f;
double 	x_l1_norm;
double 	F_B;
double 	*step; /* difference between two steps for trust region update */
double 	norm_step;
double 	xi;

double 	lambda_l2;

/* for least-square loss */
int 	saving_flag;
double 	**ATA; 		/* A^TA */
double 	*ATy; 		/* A^Ty */
double 	norm_y_sq;  /* ||y||_2^2 */
double 	*Ax;		/* Ax */
double 	*AI_x;

int 	maxback;
int 	sameOrthant;
int 	sameOrthant_prev;

Node_t *X_row;
Node_t *X_col;

char 	linesearch_flag[40];
char 	type_iteration[40];

int 	*col_ptr;
int 	*col_ptr_head;
double 	*mean_column;
double 	*min_column;
double 	*max_column;

static 	int 	max_line_size 	= 1024;
static 	char 	*line 			= NULL;
//static 	char 	*column_titles 	= " Iter      f       |x|        F      |phi|   |beta|    nnz | type    nI  TRrad    flag    its residual  target    nV  nVmax    |d|   | bak   alpha    flag    nV  | phi%% bet%% btRat spRat | btSec spSec itSec";
static 	char 	*column_titles 	= " Iter      f       |x|        F      |phi|   |beta|    nnz | type    nI  TRrad    flag  rsType   its residual  target    nV  nVmax    |d|   | bak   alpha    flag    nV  |";
static 	char 	*double_hline  	= "=====================================================================================================================";
static 	char 	*hline 			= "---------------------------------------------------------------------------------------------------------------------";
static 	char 	*hline_1 		= "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------";

/* global output FaRSA structure */
Output_FaRSA_t 	output_farsa;

/* global parameter structure for setting information */
Parameter_t 	param;

/* global problem for simplifying argument in coding */
Problem_t 		problem;


/* call farsa routine to solve target problem */
Output_FaRSA_t farsa( int argc, char **argv, Input_FaRSA_t *input_farsa );

/* Load information from parser arguments */
void 		parse_command_line( int argc, char **argv );

/* Load information from profile */
void 		load_setting_profile();

/* read line from input, from liblinear */
char* 		get_line( FILE *input );

/* initialize rest parameters */
void 		initialize_rest_param( struct Input_FaRSA input_farsa );

/* read data by the libsvm format */
void 		read_sparse_problem( struct Problem *prob, int is_train  );

/* transpose libsvm format data */
void 		transpose( struct Problem *prob );

/* set beta, phi and grad_F */
void 		setBetaPhi();

/* swap a-th entry and b-th entry in I */
void 		swap( int *a, int *b );

/* set beta, phi and grad_F */
void 		setBetaPhi();

/* select beta */
void 		select_beta();

/* select phi */
void 		select_phi();

/* quick select K elements with largest magnitude from vector */
double 		quick_select( double *vector, int K, int lo, int hi, int *indexes );

/* quick sort with duplicate keys by three-way partition */
void  		quick_sort( int* array, int lo, int hi );

int 		partition_descend( double *vector, int lo, int hi, int *indexes );

/* optimized routine for logistic cost function */
void 		logistic_loss( struct Problem *prob  );

/* update some intermediate variables of logistic loss */
void 		logistic_setExpTerm(  double *x, double *wTx, double *sigmoid, double *ysigmoid, double *expterm );

/* calculate logistic loss function value */
Func_output_t logistic_func( double *x, double *expterm );

/* set data for phi-step subproblem */
void 		logistic_setXF();

/* logistic loss hess vector product */
void 		logistic_hessVecProd( double *v, double *hv );


/* optimized routine for elastic_net */
void 		elastic_net_loss( struct Problem *prob );

/* update Ax */
void 		least_sqaure_setAx( double *x );

/* calculate least-square loss function value */
Func_output_t least_square_func( double *x );

/* least square loss hess vector product */
void 		least_square_hessVecProd( double *v, double *hv );

/* least square loss set data for phi-step subproblem */
void 		least_square_setXF();

/* reduced space CG solver */
Output_sp_t CGsolver( Input_CG_t input_cg );

/* reduced space CD solver */
Output_sp_t CDsolver_DL_diffd_logloss( Input_CD_t input_cg );

/* Projected gradient descent */
Proj_output_t project( double *x, double *d, double alpha, double *proj_x );

/* calculate sparsity in solution */
double 		calculate_sparsity( double *x, int n );

/* Sparse operators start */
/* Operator 1: the inner product of one vector and a sparse vector with the same dimension */
double 		dot_ds( const double *v, const Node_t *x );

/* Operator 2: the inner product of two sparse vectors */
double 		dot_ss( const Node_t *v_1, const Node_t *v_2 );

/* n5 operators start */
/* Operator 7: n5 type inner product */
double 		dot_n5( double *v1, double *v2, int n );

/* Operator 8: n5 type l1 norm */
double 		l1_n5( double *v, int n );


/* Output function based on print level*/
void 		output();

/* Output header */
void 		output_header();

/* Output problem data */
void 		output_problem_data();

/* Output parameter controls */
void 		output_param();

/* output final x */
void 		output_final_x();

/* test routine */
Test_output_t test( Problem_t *problem_test );

/* predict label of one instance by logisitc regression */
double 		logistic_predict_one_instance( double *x, const Node_t *d,double bias );

/* output test results */
void 		output_test( Test_output_t test_result );

/* print tools */
void 		print( double *vector, int n );

void 		printe( double *vector, int n );

void 		printInt( int *vector, int n );

void 		print_sparse_vec( const Node_t *start_node );

/* copy the elements in vector2 to vector1 */
void 		copy( double *vector1, double *vector2, int n );

#endif