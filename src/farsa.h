/* Author: TC, FEC, DPR */
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
#define LEAST_SQUARE	2

/* macro for sub-problem sover */
#define USE_ADAPTIVE 	0
#define USE_DIRECT 		1
#define USE_CG 			2
#define USE_CD_QP 		3

#define MAX_VALUE 		1E7

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define MAX( x, y ) ( ((x) > (y)) ? (x) : (y) )
#define MIN( x, y ) ( ((x) < (y)) ? (x) : (y) )
#ifndef FALSE
#define FALSE 			0
#endif
#ifndef TRUE
#define TRUE 			1
#endif


/* Sparse data node that contains the row index, column index and the value of data point*/
typedef struct Node_rc {

	int 	r_index; /* r_index is the index of row that is index of sample */
	int 	c_index; /* c_index is the index of column that is index of feature */
	double 	value;

} Node_rc_t;


/* Problem structure to save information */
typedef struct Problem {

	int 		num_features;
	int 		num_samples;
	double 		*y;
	Node_rc_t 	**rX; 
	Node_rc_t 	**cX;
	Node_rc_t 	**cXS;
	Node_rc_t 	**rXS;

} Problem_t;


/* Parameter structure */
typedef struct Parameter {

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

	char 	prob_name[1024];
	char 	profile_file[1024];
	char 	data_format[1024];
	char 	output_file[1024];
	char 	name[1024];

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


/* Input argument for farsa */
typedef struct Input_FaRSA{

	int 	n;
	double 	( *func )( double * );
	double 	*( *grad_f )( double * );
	double 	*( *hessVec )( double * );

} Input_FaRSA_t;


typedef struct Output_FaRSA{

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
typedef struct Proj_output{

	double 	*project_vector;	
	int 	same_sign;
	int 	nV;

} Proj_output_t;


/* Function value result */
typedef struct Func_output{
	
	double 	f_value;
	double	F_value;
	double 	x_l1_norm;

} Func_output_t;

/* Subspace solver part starts */

/* structure of basic dog-leg search sub-problem solver */
typedef struct Input_direct {

	int 	*I;
	int 	*I_p;
	int 	*I_n;
	int 	nI;
	int 	nVmaxAllow;
	int 	max_direct_iter;
	int 	n;
	double 	*grad_F;
	double 	eta_r;
	double 	rsType;
	double 	fracViol;
	double 	TRradius;
	double 	HvProds;
	double 	regularize;
	double 	*x;
	double 	*( *hessVecProd )( double *);

} Input_direct_t;

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
	double 	eta_r;
	double 	rsType;
	double 	fracViol;
	double 	TRradius;
	double 	HvProds;
	double 	*( *hessVecProd )( double *);

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


typedef struct Output_CD{

	int 	nV;
	int 	iter;
	double 	*d_full;
	double 	dir_der;
	char 	sub_prob_flag[1024];

} Output_CD_t;

/* Subspace solver part ends */

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
double 	dm;

/* saving the diagonal entries of D in X^TDX for logistic loss */
double 	*diagonal;
double 	*y;
double 	*y_B;
double 	*d;
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

int 	maxback;
int 	sameOrthant;
int 	sameOrthant_prev;

Node_rc_t *X;
Node_rc_t *X_col;
Node_rc_t *X_row_S;


char 	linesearch_flag[40];
char 	type_iteration[40];

int 	*col_ptr;
int 	*col_ptr_head;
double 	*mean_column;
double 	*min_column;
double 	*max_column;

static 	int 	max_line_size 	= 0;
static 	char 	*line 			= NULL;
//static 	char 	*column_titles 	= " Iter      f       |x|        F      |phi|   |beta|    nnz | type    nI  TRrad    flag    its residual  target    nV  nVmax    |d|   | bak   alpha    flag    nV  | phi%% bet%% btRat spRat | btSec spSec itSec";
static 	char 	*column_titles 	= " Iter      f       |x|        F      |phi|   |beta|    nnz | type    nI  TRrad    flag  rsType   its residual  target    nV  nVmax    |d|   | bak   alpha    flag    nV  |";
static 	char 	*double_hline  	= "=====================================================================================================================";
static 	char 	*hline 			= "---------------------------------------------------------------------------------------------------------------------";
static 	char 	*hline_1 		= "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------";


/* global problem for simplifying argument in coding */
Problem_t 	problem;

/* global parameter structure for setting information */
Parameter_t param;

/* global output FaRSA structure */
Output_FaRSA_t output_farsa;

/* Load information from parser arguments */
void 		parse_command_line( int argc, char **argv );

/* Load information from profile */
void 		load_setting_profile();

char* 		get_line( FILE *input );

/* select the correct way to read problem */
void 		read_problem( struct Problem *prob  );

/* read data by the libsvm format */
void 		read_sparse_problem( struct Problem *prob  );

/* transpose libsvm format data */
void 		transpose( struct Problem *prob );

/* read dense type data */
void 		read_dense_problem( struct Problem *prob );

/* read libsvm format data with scaling techniques */
void 		read_scaling_problem( struct Problem *prob  );

/* call farsa routine to solve target problem */
Output_FaRSA_t farsa( int argc, char **argv, Input_FaRSA_t *input_farsa );

/* initialize rest parameters */
void 		initialize_rest_param( struct Input_FaRSA input_farsa );

/* reduced space sub-problem solver */
Output_sp_t Directsolver( Input_direct_t input_direct );

/* reduced space CG solver */
Output_sp_t CGsolver( Input_CG_t input_cg );

/* reduced space CD solver */
Output_sp_t CDsolver_DL_diffd_logloss( Input_CD_t input_cg );

/* Projected gradient descent */
Proj_output_t project( double *x, double *d, double alpha );

/* swap a-th entry and b-th entry in I */
void 		swap( int *a, int *b );

/* set beta, phi and grad_f */
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

void 		logistic_loss( struct Problem *prob  );

/* calculate logistic loss function value */
Func_output_t logistic_func( double *x, double *expterm );

/* update some intermediate variables of logistic loss */
void 		logistic_setExpTerm(  double *x, double *wTx, double *sigmoid, double *ysigmoid, double *expterm );

/* logistic loss hess vector product */
double 		*logistic_hessVecProd( double *v );

void 		logistic_setXF();

/* generic objective function routine */
void 		generic_loss( struct Input_FaRSA *input_farsa );
	
/* calculate sparsity in solution */
double 		calculate_sparsity( double *x, int n );

/* Sparse operators start */
/* Operator 1: the inner product of one vector and a sparse vector with the same dimension via row */
double 		dot_r( const double *v, const Node_rc_t *x );

/*  Operator 2: the inner product of one vector and a sparse vector with the same dimension via column */
double 		dot_c( const double *v, const Node_rc_t *x );


double 		dot_r_divide_num( const double *v, const Node_rc_t *x, const double num  );
 
/*  Operator 3: l2 norm square by row */
double 		norm2_sq_r( const Node_rc_t *x );

/*  Operator 4: l2 norm square by column */
double 		norm2_sq_c( const Node_rc_t *x );

/*  Operator 5: inner product of v and subvector of the sparse vector given index S by row. */
double 		dot_S_r( const double *v, const Node_rc_t *x, int *S );

/*  Operator 6: inner product of v and subvector of the sparse vector given index S by column. */
double 		dot_S_c( const double *v, const Node_rc_t *x, int *S );

/* Sparse operators end */


/* n5 operators start */
/* Operator 7: n5 type inner product */
double 		dot_n5( double *v1, double *v2, int n );

/* Operator 8: n5 type l1 norm */
double 		l1_n5( double *v, int n );

/* n5 operators end */

/* Free everything to avoid crash */
void 		free_all();

/* Output function based on print level*/
void 		output();

/* Output header */
void 		output_header();

/* Output problem data */
void 		output_problem_data();

/* Output parameter controls */
void 		output_param();

/* print tools */
void 		print( double *vector, int n );

void 		printe( double *vector, int n );

void 		printInt( int *vector, int n );

/* copy the elements in vector2 to vector1 */
void 		copy( double *vector1, double *vector2, int n );



#endif 