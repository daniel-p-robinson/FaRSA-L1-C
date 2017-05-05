/* Authors: TC, DPR, FEC */
#include "farsa.h"

Output_FaRSA_t farsa( int argc, char **argv, Input_FaRSA_t *input_farsa ){

	clock_t begin;
	clock_t end;
	double run_time;
	/* get arguments from parser */
	parse_command_line( argc, argv );

	/* load setting from profile */
	load_setting_profile();

	// fprintf( stdout, "%d\n", param.sub_solver );

	/* initialize rest parameters */
	/* load dataset is contained here. */
	initialize_rest_param( *input_farsa );

	// fprintf( stdout, "Load data complete\n" );

	/* if loss function has optimized version, then directly call optimized routine
	  if not, run the generic routine. */
	begin = clock();

	switch ( param.objective_function_type ) {

		case GENERIC_FUNC:
		
			if( param.print_level >= 1 && param.print_level != 3 && param.print_level != 4 ){
				fprintf( stdout, "Personalized loss function...\n" );
			}
			generic_loss( input_farsa );
			break;
		
		case LOGISTIC_LOSS:
			
			if( param.print_level >= 1 && param.print_level != 3 && param.print_level != 4 ){
				fprintf( stdout, "Logistic loss plus l1 regularizer...\n" );
			}
			logistic_loss( &problem );
			break;
		
		case LEAST_SQUARE:
		
			if( param.print_level >= 1 && param.print_level != 3 && param.print_level != 4 ){
				fprintf( stdout, "Least Square Loss plus l1 regularizer...\n" );
			}
			break;

		default:

			if( param.print_level >= 1 ){
				fprintf( stdout, "objective function type is invalid \n" );
			}
			exit(0);	

	}

	end 					= clock();
	run_time 				= (double)(end - begin) / CLOCKS_PER_SEC;
	output_farsa.run_time 	= run_time;
	/* output based on print level */
	output();

	return output_farsa;
}

/****************** functions that can be used for any function including special optimized loss function ***/
Proj_output_t project( double *x, double *d, double alpha ){

	int 			i;
	int 			same_sign;
	int 			nV;
	double 			*x_proj_line;
	Proj_output_t 	proj_output;
	
	same_sign 		= TRUE;
	x_proj_line 	= Malloc( double, n );
	nV 				= 0;

	for( i = 0; i < n; i++ ){
	
		x_proj_line[i] = x[i] + alpha * d[i];
	
		if( I_p[i] == 1 && x_proj_line[i] < 0 ){
	
			x_proj_line[i] = 0.0;
			same_sign = FALSE;
			nV++;

		}else if( I_n[i] == 1 && x_proj_line[i] > 0 ){
	
			x_proj_line[i] = 0.0;
			same_sign = FALSE;
			nV++;

		}

	}
	
	proj_output.project_vector 	= x_proj_line;
	proj_output.same_sign 		= same_sign;
	proj_output.nV 				= nV;

	return proj_output;
}

void setBetaPhi() {	

	int 			i;
	int 			*p_phi;
	int 			*p_beta;

	p_phi 			= I_phi;
	p_beta 			= I_beta;

	norm_beta 		= 0.0;
	norm_phi 		= 0.0;

	beta 			= Malloc( double, n );
	phi 			= Malloc( double, n );

	nnz 			= 0;

	for( i = 0; i < n; i++ ){
		beta[i] = 0.0;
		phi[i] = 0.0;
	}


	for ( i = 0; i < n; i++ ){

		if( x[i] == 0 ){
			I_z[i] = 1;
			I_p[i] = 0;
			I_n[i] = 0;
			if ( grad_f[i] < -lambda ){
				beta[i] = grad_f[i] + lambda;			
			}else if( grad_f[i] > lambda ){
				beta[i] = grad_f[i] - lambda;
			}
			/* use l2 norm at first */
			norm_beta += beta[i] * beta[i];
			phi[i] = 0.0;

			if ( beta[i] != 0.0 ) {
				nI_beta++;
				*p_beta = i;
				p_beta++;
			}

		}else if( x[i] > 0 ){
			phi[i] = grad_f[i] + lambda;
			if( phi[i] > 0 ){
				phi[i] = MIN( phi[i], MAX( x[i], grad_f[i] - lambda ) );
			}
			norm_phi += phi[i] * phi[i];
			I_p[i] = 1;
			I_z[i] = 0;
			I_n[i] = 0;
			grad_F[i] = grad_f[i] + lambda; 
			if( phi[i] != 0.0 ){
				nI_phi++;
				*p_phi = i;
				p_phi++;
			}
			
			nnz++;

		}else{
			phi[i] = grad_f[i] - lambda;
			if( phi[i] < 0 ){
				phi[i] = MAX( phi[i], MIN( x[i], grad_f[i] + lambda ) );
			}
			norm_phi += phi[i] * phi[i];
			I_n[i] = 1;
			I_z[i] = 0;
			I_p[i] = 0;
			grad_F[i] = grad_f[i] - lambda;
			if( phi[i] != 0.0 ){
				nI_phi++;
				*p_phi = i;
				p_phi++; 
			}

			nnz++;
		}

	}

	*p_phi 	= -1;
	p_phi 	= I_phi;
	*p_beta = -1;
	p_beta	= I_beta;

	norm_phi 	= sqrt( norm_phi );
	norm_beta 	= sqrt( norm_beta );
	// fprintf(stdout, "norm of grad_f : %e\n", sqrt( dot_n5( grad_f, grad_f, n ) ) ) ;
	// fprintf(stdout, "nI_beta: %d\n", nI_beta );

}

/* select beta */
void select_beta() {
	if ( param.betaFrac >= 1.0 ) {
		Beta 		= beta;
		norm_Beta 	= norm_beta;
		return ;
	} else {	

		int i;
		
		beta_I 				= Malloc( double, nI_beta );
		I_beta_selected 	= Malloc( int, nI_beta );

		for ( i = 0; i < nI_beta; i++ ) {
			beta_I[i] 			= beta[ I_beta[i] ];
			I_beta_selected[i] 	= I_beta[i];
		}

		nI_beta_selected 	= MIN( nI_beta, MAX( 1, floor( (double) n * param.betaFrac ) ) );
		quick_select( beta_I, nI_beta_selected - 1, 0, nI_beta - 1, I_beta_selected );
		quick_sort( I_beta_selected, 0, nI_beta_selected - 1 );

		/* update Beta, norm_Beta */
		for ( i = 0; i < n; i++ ) {
			Beta[i] 	= 0.0;
		}

		norm_Beta 	= 0.0;
		for ( i = 0; i < nI_beta_selected; i++ ) {
			Beta[I_beta_selected[i]]	= beta[I_beta_selected[i]];
			norm_Beta += Beta[I_beta_selected[i]] * Beta[I_beta_selected[i]];
		}

		norm_Beta 	= sqrt( norm_Beta );

	}
}

/* select phi */
void select_phi() {

	if ( param.phiFrac >= 1.0 ) {
		I_phi_selected 		= I_phi;
		nI_phi_selected 	= nI_phi;
		return ;

	} else {

		int i;
		phi_I 			= Malloc( double, nI_phi );
		I_phi_selected 	= Malloc( int, nI_phi );

		for ( i = 0; i < nI_phi; i++ ) {
			phi_I[i] 			= phi[I_phi[i]];
			I_phi_selected[i]	= I_phi[i];
		}
		nI_phi_selected = MIN( nnz, MAX( ceil( (double) nnz * param.phiFrac ), 5000 ) );
		quick_select( phi_I, nI_phi_selected - 1, 0, nI_phi - 1, I_phi_selected );
		quick_sort( I_phi_selected, 0, nI_phi_selected - 1 );
	
	}

}

/* quick select largest K from vector */
double quick_select( double *vector, int K, int lo, int hi, int *indexes ) {

	int i;
	while ( lo <= hi ) {
		i = partition_descend( vector, lo, hi, indexes );
		if ( i < K ) {
			lo = i + 1;
		} else if ( i > K ) {
			hi = i - 1;
		} else {
			return vector[K];
		}
	}
	return vector[K];

}

int partition_descend( double *vector, int lo, int hi, int *indexes ) {
	int 	i, j;
	double 	pivot;
	double 	temp;
	int  	temp_idx;

	pivot 	= fabs(vector[lo]);
	i 		= lo;
	j 		= hi;

	while( TRUE ) {

		while ( fabs(vector[i]) >= pivot && i <= hi ) {
			i++;
		}

		while ( fabs(vector[j]) < pivot ) {
			j--;
		}

		if ( i >= j ) {
			break;
		}

		temp 		= vector[j];
		vector[j] 	= vector[i];
		vector[i]	= temp;
		temp_idx 	= indexes[j];
		indexes[j]	= indexes[i];
		indexes[i] 	= temp_idx;
	}

	temp 		= vector[j];
	vector[j]	= vector[lo];
	vector[lo]	= temp;
	temp_idx 	= indexes[j];
	indexes[j]	= indexes[lo];
	indexes[lo] 	= temp_idx;
	
	return j;

}

/* quick sort with duplicate keys by three-way partition */
void quick_sort( int* array, int lo, int hi ) {
	
	int median_idx;
	int pivot_left_idx;
	int pivot_right_idx;
	int pivot;
	int i;
	if ( hi <= lo ) {
		return;
	}
	median_idx 		= lo + ( hi - lo ) / 2;
	pivot_left_idx 	= lo;
	pivot_right_idx = hi;
	swap( &array[lo], &array[median_idx] );

	pivot 			= array[lo];
	i 				= lo;

	while ( i <= pivot_right_idx ) {
		if ( array[i] < pivot ) {
			swap( &array[pivot_left_idx], &array[i] );
			pivot_left_idx++;
			i++;
		} else if ( array[i] > pivot ) {
			swap( &array[i], &array[pivot_right_idx] );
			pivot_right_idx--;
		} else {
			i++;
		}		
	}

	quick_sort( array, lo, pivot_left_idx - 1 );
	quick_sort( array, pivot_right_idx + 1, hi );

}

/*************************************************************/

/***********  Personalized Loss Section Start ****************/

void generic_loss ( struct Input_FaRSA *input_farsa ) {

	int 			i, j;
	int 			iter;
	double 			initial_F;
	double 			suf_descent_1;

	Proj_output_t 	proj_output;
	Input_CG_t		input_cg;
	Output_sp_t 	output_cg;


	iter 			= 0;
	F 				= input_farsa->func(x);	
	initial_F 		= F;
	grad_f 			= input_farsa->grad_f( x );

	grad_eval++;
	f_eval++;



	while( TRUE ) {

		setBetaPhi();

		if( iter == 0 ){
			
			norm_phi0 = norm_phi;
			norm_beta0 = norm_beta;

			if( param.termType == 2 ){
				ttol = MAX( tol_absolute, tol_relative * MAX( norm_beta0, norm_phi0 ) );
			}else if( param.termType == 1 ){
				ttol = tol_absolute;
			}

			output_farsa.x_initial 			= x;
			output_farsa.norm_beta_initial 	= norm_beta0;
			output_farsa.norm_phi_initial 	= norm_phi0;
			output_farsa.ttol 				= ttol;
		
		}

		// termination
		if( MAX( norm_beta, norm_phi ) <= ttol ){
			/* set output structure */
			output_farsa.iter 				= iter;
			output_farsa.beta_iter 			= beta_iter;
			output_farsa.phi_iter 			= phi_iter;
			output_farsa.x_final 			= x;
			output_farsa.norm_beta_final 	= norm_beta;
			output_farsa.norm_phi_final 	= norm_phi;
			output_farsa.F_final 			= F;
			output_farsa.error 				= MAX( norm_beta, norm_phi );
			output_farsa.f_eval 			= f_eval;
			output_farsa.grad_eval 			= grad_eval;
			output_farsa.hessVec_eval 		= hessVec_eval;
			output_farsa.term 				= 1;
			output_farsa.zero_perc 			= calculate_sparsity( x, n );
			break;
		}

		/* if the number of iteration reached the maximum iteration, break */
		if( iter >= param.max_iter ){

			/* set output structure */
			output_farsa.iter 				= param.max_iter;
			output_farsa.beta_iter 			= beta_iter;
			output_farsa.phi_iter 			= phi_iter;
			output_farsa.x_final 			= x;
			output_farsa.norm_beta_final 	= norm_beta;
			output_farsa.norm_phi_final 	= norm_phi;
			output_farsa.F_final 			= F;
			output_farsa.error 				= MAX( norm_beta, norm_phi );
			output_farsa.f_eval 			= f_eval;
			output_farsa.grad_eval 			= grad_eval;
			output_farsa.hessVec_eval 		= hessVec_eval;
			output_farsa.term 				= 0;
			output_farsa.zero_perc 			= calculate_sparsity( x, n );
			break;

		}
		
		if( norm_beta <= Gamma * norm_phi ){


			
			phi_iter++;

			iter_type = 1;
			
			input_cg.n 			= n;
			input_cg.nI 		= nI_phi;
			input_cg.x 			= x;
			input_cg.I 			= I_phi;
			input_cg.I_p 		= I_p;
			input_cg.I_n 		= I_n;
			input_cg.eta_r 		= param.eta_r;
			input_cg.grad_F 	= grad_F;
			input_cg.nVmaxAllow = ceil( MAX((double)param.nVmaxAllow,input_cg.fracViol*(double)nI_phi) );
			input_cg.TRradius 	= TRradiusPhi;
			input_cg.maxCG_iter = param.max_CG_iter;
			input_cg.fracViol 	= param.fracViol;
			input_cg.rsType 	= rsType;
			input_cg.hessVecProd = input_farsa->hessVec;

			output_cg = CGsolver( input_cg );

			d = output_cg.d_full;

			if( output_cg.nV == 0 && norm_beta < tol_absolute ){
				rsType = 2;
			}else{
				rsType = 1;
			}

			/* Perform a line search */
			j = 0;
			alpha = 1.0;
			dirDer = output_cg.dir_der;
			suf_descent_1 = param.eta * dirDer; 
			F_old = F;
			proj_output = project( x, d, alpha );
			x_linesearch = proj_output.project_vector;
			sameOrthant = proj_output.same_sign;
			F = input_farsa->func(x_linesearch);
			f_eval++;

			while(1){
				
				if( sameOrthant && F - F_old <=  suf_descent_1 * alpha ){
					
					/* set x, and grad_f */
					for( i = 0; i < n; i++ ){
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						x[i] = x_linesearch[i];
					}
					norm_step = sqrt( norm_step );
					break;
				}

				if( !sameOrthant && F < F_old ){
					// set x, and grad_f
					for( i = 0; i < n; i++ ){
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						x[i] = x_linesearch[i];
					}
					norm_step = sqrt( norm_step );
					break;
				}

				if ( j > maxback ){
					// set x, and grad_f
					for( i = 0; i < n; i++ ){
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						x[i] = x_linesearch[i];
					}
					norm_step = sqrt( norm_step );
					break;
				}

				alpha *= xi;
				proj_output = project( x, d, alpha );
				x_linesearch = proj_output.project_vector;
				sameOrthant = proj_output.same_sign;
				F = input_farsa->func(x_linesearch);
				f_eval++;
				j++;

			}	
			
			grad_f = input_farsa->grad_f( x_linesearch );	
			grad_eval++;


		}else{

			beta_iter++;
			iter_type = 0;
			
			select_beta();

			double TRbeta_norm_beta = -( TRradiusBeta / norm_Beta );

			for ( i = 0; i < n; i++ ){
				d[i] = TRbeta_norm_beta * Beta[i];
				x_linesearch[i] = x[i] + d[i];
				if( d[i] > 0 ){
					grad_F[i] = grad_f[i] + lambda;
				}else if( d[i] < 0 ){
					grad_F[i] = grad_f[i] - lambda;
				}else{
					grad_F[i] = 0.0;
				}
			}
			grad_eval++;

			double normd = TRradiusBeta;
			/* Perform line search to get updated x */
			j = 0;
			alpha = 1.0;
			dirDer = dot_n5( grad_F, d, n );

			F_old = F;

			F = input_farsa->func(x_linesearch);
			f_eval++;

			double suf_descent_1 = param.eta * dirDer; 
			double tmp = 0.0; // used for update grad_f

			/* perform line search */
			while( 1 ){
				if( F - F_old <= alpha * suf_descent_1 ){
					/* set x, and grad_f */
					for( i = 0; i < n; i++ ){
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						x[i] = x_linesearch[i];
						// grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
					}
					norm_step = sqrt( norm_step );
					break;
				}
				if ( j > maxback ){
					/* set x, and grad_f */
					for( i = 0; i < n; i++ ){
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						x[i] = x_linesearch[i];
						// grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
					}
					norm_step = sqrt( norm_step );
					break;
				}
				j++;
				alpha *= xi;	

				/* set new trial step */
				for( i = 0 ; i < n; i++ ){
					x_linesearch[i] = x[i] + alpha * d[i];
				}

				F = input_farsa->func(x_linesearch);
				f_eval++;
			}

			grad_f = input_farsa->grad_f( x_linesearch );
			grad_eval++;

			// fprintf(stdout, "%f\n", F);

		}

		d 		= Malloc( double, n );
		nI_phi 	= 0;
		nI_beta	= 0;

		/* Set trust-region radius for next iteration. */
		if( iter_type == 1 ){
			TRradiusPhi = MAX( 1E-3, MIN( 1E3, 10 * norm_step ) );
		}else{
			TRradiusBeta = MAX( 1E-5, MIN( 1.0, norm_step ) );
		}
		norm_step = 0.0;		
		iter++;
	
	}


}


/***********  Personalized Loss Section End   ****************/

/***********  Logistic Regression Loss Section Start *********/
/**
 *	The routine contains 
 */
void logistic_loss( Problem_t *prob ){

	int 			i;
	int 			j;
	int 			nDifforthant;
	int 			nV;
	int 			nVmax;

	double 			*sigmoid;
	double 			*ysigmoid;
	double 			*expterm;
	double 			*wTx;
	double 			suf_descent_1;
	double 			initial_F;
	double 			norm_d;
	double 			TRbeta_norm_beta;
	double 			tmp;
	double 			sp_residual;
	double 			sp_target_residual;
	FILE 			*fo;
	char 			linesearch_flag[40];

	Func_output_t	func_output;
	Proj_output_t 	proj_output;
	Input_CG_t 		input_cg;
	Input_direct_t	input_dt;
	Input_CD_t 		input_cd;
	Output_sp_t 	output_cp;

	/* help_vector for logistic loss is sigmoid vector */
	diagonal 		= Malloc( double, m );
	sigmoid 		= Malloc( double, m );
	ysigmoid 		= Malloc( double, m );
	expterm 		= Malloc( double, m );
	wTx 			= Malloc( double, m );

	/* Set labels and Initialize some expterm ysigmoid */
	for ( i = 0; i < m; i++ ) {
		if ( problem.y[i] > 0 ) {
			y[i] 	= 1.0;
		} else {
			y[i] 	= -1.0;
		}
		expterm[i] 	= 1.0;
		sigmoid[i] 	= 0.5;
		ysigmoid[i] = y[i] * sigmoid[i];
		diagonal[i] = 0.25;
	}


	/* Continue initialize */
	for( i = 0 ; i < n; i++ ){
		x[i] 		= 0.0;
		grad_f[i] 	= -1.0 / dm * dot_r( ysigmoid, prob->cX[i] );
	}
	grad_eval++;

	func_output 	= logistic_func( x, expterm );
	F 				= func_output.F_value;
	f_eval++;


	initial_F 				= F;
	output_farsa.F_initial 	= initial_F;
	output_farsa.lambda 	= lambda;

	/* Initially set beta phi */
	setBetaPhi();

	norm_phi0 		= norm_phi;
	norm_beta0 		= norm_beta;

	/* set target termination tolerance */
	if( param.termType == 2 ){
		ttol = MAX( tol_absolute, tol_relative * MAX( norm_beta0, norm_phi0 ) );
	}else if( param.termType == 1 ){
		ttol = tol_absolute;
	}

	output_farsa.x_initial 			= x;
	output_farsa.norm_beta_initial 	= norm_beta0;
	output_farsa.norm_phi_initial 	= norm_phi0;
	output_farsa.ttol 				= ttol;
	
	/* print to output file if print level = 4 */
	if ( param.print_level == 4 ) {

		fo = fopen( param.output_file, "w" );
		if ( fo == NULL ) {
			fprintf( stdout, "Can not open %s\n", param.output_file );
			exit(0);
		}

		fprintf( fo, "%s\n", hline );
		fprintf( fo, "|                  A Fast Reduced Space Algorithm (FaRSA) for Sparse Logistic Regression v. 1.0                     |\n" );
		fprintf( fo, "%s\n", double_hline );

		fclose( fo );

		output_problem_data();
		
		output_param();
	
	}


	/* Main loop for FaRSA for logistic loss */
	while( TRUE ) {	

		/* print based on print level */
		if ( param.print_level >= 4 ) {

			fo = fopen( param.output_file, "ab" );
			
			if ( iter % param.print_every == 0 ) {

				fprintf( fo, "%s\n", hline_1 );
				fprintf( fo, "%s\n", column_titles );
				fprintf( fo, "%s\n", hline_1 );			
			}
			
			fprintf( fo, "%5d %9.2e %8.2e %9.2e %8.2e %8.2e %5d |", iter, func_output.f_value, func_output.x_l1_norm, func_output.F_value, norm_phi, norm_beta, nnz );
			
			fclose( fo );

		}

		/* termination condition */
		/* terminate if error is smaller than the tolerance */
		if( MAX( norm_beta, norm_phi ) <= ttol ){

			output_farsa.iter 			 = iter;
			output_farsa.beta_iter 		 = beta_iter;
			output_farsa.phi_iter 		 = phi_iter;
			output_farsa.x_final 		 = x;
			output_farsa.norm_beta_final = norm_beta;
			output_farsa.norm_phi_final  = norm_phi;
			output_farsa.F_final 		 = F;
			output_farsa.f_final 		 = f;
			output_farsa.x_l1_norm_final = x_l1_norm;
			output_farsa.error 			 = MAX( norm_beta, norm_phi );
			output_farsa.f_eval 		 = f_eval;
			output_farsa.grad_eval 		 = grad_eval;
			output_farsa.hessVec_eval 	 = hessVec_eval;
			output_farsa.term 			 = 1;
			output_farsa.ttol 			 = ttol;
			output_farsa.zero_perc 		 = calculate_sparsity( x, n );
			output_farsa.status 		 = 0;
			break;

		}

		/* terminate if current iteration exceed maximum allowed iteration */
		if( iter >= param.max_iter ){
			
			output_farsa.iter 			 = param.max_iter;
			output_farsa.beta_iter 		 = beta_iter;
			output_farsa.phi_iter 		 = phi_iter;
			output_farsa.x_final 		 = x;
			output_farsa.norm_beta_final = norm_beta;
			output_farsa.norm_phi_final  = norm_phi;
			output_farsa.F_final 		 = F;
			output_farsa.f_final 		 = f;
			output_farsa.x_l1_norm_final = x_l1_norm;
			output_farsa.error 			 = MAX( norm_beta, norm_phi );
			output_farsa.f_eval 		 = f_eval;
			output_farsa.grad_eval 		 = grad_eval;
			output_farsa.hessVec_eval 	 = hessVec_eval;
			output_farsa.term 			 = 0;
			output_farsa.ttol 			 = ttol;
			output_farsa.zero_perc		 = calculate_sparsity( x, n );
			output_farsa.status 		 = 1;
			break;
		
		}

		/* if norm of beta is less than the norm of phi, then do phi step */
		if( norm_beta <= Gamma * norm_phi ){

			sprintf( type_iteration, "%s", "phi" );
			
			phi_iter++;
			
			iter_type = 1;

			select_phi();


			/* solve sub-problem based on corresponding type */
			switch ( param.sub_solver ) {

				case USE_ADAPTIVE:
				
					break;
		
				case USE_CG:

					input_cg.n 			 = n;
					input_cg.nI 		 = nI_phi_selected;
					input_cg.x 			 = x;
					input_cg.I 			 = I_phi_selected;
					input_cg.I_p 		 = I_p;
					input_cg.I_n 		 = I_n;
					input_cg.eta_r 		 = param.eta_r;
					input_cg.grad_F 	 = grad_F;
					nVmax 				 = ceil( MIN((double)param.nVmaxAllow,param.fracViol*(double)nI_phi_selected) );
					input_cg.nVmaxAllow  = nVmax;
					input_cg.TRradius 	 = TRradiusPhi;
					input_cg.maxCG_iter  = param.max_CG_iter;
					input_cg.fracViol 	 = param.fracViol;
					input_cg.rsType 	 = rsType;
					input_cg.iter 		 = iter;
					input_cg.hessVecProd = &logistic_hessVecProd;
					
					logistic_setXF();

					output_cp 			 = CGsolver( input_cg );


					free(problem.cXS);
					free(problem.rXS);
					free( X_row_S );
					problem.cXS = NULL;
					problem.rXS = NULL;
			
					break;
				
				/* Compute the search direction using CD on the quadratic model. */
				case USE_CD_QP:
					
					input_cd.n 			 = n;
					input_cd.m 			 = m;
					input_cd.nI 		 = nI_phi_selected;
					input_cd.I 			 = I_phi_selected;
					input_cd.grad_F 	 = grad_F;
					input_cd.eta_r 		 = param.eta_r;
					nVmax 				 = ceil( MIN((double)param.nVmaxAllow,param.fracViol*(double)nI_phi_selected) );
					input_cd.nVmaxAllow  = nVmax;
					input_cd.TRradius 	 = TRradiusPhi;		
					input_cd.help_vector = diagonal;			
					input_cd.rsType 	 = rsType;
					input_cd.regularize  = 1E-8;
					input_cd.grad_F_I 	 = Malloc( double, nI_phi_selected );
					input_cd.x_I 		 = Malloc( double, nI_phi_selected );
					input_cd.x_I_pos 	 = Malloc( int, nI_phi_selected );
					input_cd.x_I_neg 	 = Malloc( int, nI_phi_selected );
					input_cd.hessVecProd = &logistic_hessVecProd;
					input_cd.fracViol 	 = param.fracViol;
					input_cd.max_CD_iter = param.max_CD_iter;

 					for ( i = 0; i < nI_phi_selected; i++ ){
						input_cd.grad_F_I[i] = grad_F[ I_phi_selected[i] ];
						input_cd.x_I[i] 	 = x[ I_phi_selected[i] ];
						if ( input_cd.x_I[i] > 0 ) {
							input_cd.x_I_pos[i] = 1;
						} else if ( input_cd.x_I[i] < 0 ) {
							input_cd.x_I_neg[i] = 1;
						}
					}

					logistic_setXF();

					output_cp 			 = CDsolver_DL_diffd_logloss( input_cd );

					free(problem.cXS);
					free(problem.rXS);
					free( X_row_S );
					problem.cXS = NULL;
					problem.rXS = NULL;

					break;


				default:

					break;	

			}


			/* Perform a line search */
			nV 					= output_cp.nV;
			j 					= 0;
			alpha 				= 1.0;
			nDifforthant 		= 0;
			dirDer 				= output_cp.dir_der;
			d 					= output_cp.d_full;
			sp_residual 		= output_cp.res;
			sp_target_residual 	= output_cp.res_target;
			norm_d 				= output_cp.norm_d;
			nVmax 				= output_cp.nVmax;
			suf_descent_1 		= param.eta * dirDer; 
			F_old 				= F;
			
			proj_output 		= project( x, d, alpha );
			nDifforthant 		= proj_output.nV;
			x_linesearch 		= proj_output.project_vector;
			sameOrthant 		= proj_output.same_sign;
			sameOrthant_prev	= TRUE;
			logistic_setExpTerm( x_linesearch, wTx, sigmoid, ysigmoid, expterm );
			func_output 		= logistic_func( x_linesearch, expterm );
			F 					= func_output.F_value;
			f 					= func_output.f_value;
			x_l1_norm 			= func_output.x_l1_norm;
			sub_iters 			= output_cp.iter;
			sprintf( linesearch_flag, "%s", output_cp.sub_prob_flag );
			f_eval++;

			if( nV == 0 && norm_beta < tol_absolute ){
				rsType = 2;
			}else{
				rsType = 1;
			}

			/* while loop for line search */
			while( TRUE ) {
				
				if( !sameOrthant && F < F_old ){
					/* set x, and grad_f */
					for( i = 0; i < n; i++ ){
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						x[i] = x_linesearch[i];
						grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
					}
					grad_eval++;
					norm_step = sqrt( norm_step );
					sprintf( linesearch_flag, "%s", "orthant" );
					break;
				}

				if( sameOrthant && F - F_old <=  suf_descent_1 * alpha ){
					
					if ( sameOrthant_prev || alpha > MIN( 1E-4, ttol ) ) {
						
						for( i = 0; i < n; i++ ){
							step[i] = x_linesearch[i] - x[i];
							norm_step += step[i] * step[i];
							x[i] = x_linesearch[i];
							grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
						}
					
						grad_eval++;
						norm_step = sqrt( norm_step );
						sprintf( linesearch_flag, "%s", "descent" );
						break;
					
					} else {

						alpha_B 		= MAX_VALUE;
						y_B 			= Malloc( double, n );

						for ( i = 0; i < n; i++ ) {

							/* x[i] / d[i] < 0.0 */
							if ( d[i] < 0.0 && x[i] > 0.0 ) {
								alpha_B = MIN( alpha_B, -x[i] / d[i] );
							}

							if ( d[i] > 0.0 && x[i] < 0.0 ) {
								alpha_B = MIN( alpha_B, -x[i] / d[i] );
							}

						}


						for ( i = 0; i < n; i++ ) {							
							y_B[i] 		= x[i] + alpha_B * d[i];
						}

						logistic_setExpTerm( y_B, wTx, sigmoid, ysigmoid, expterm );
						func_output 		= logistic_func( y_B, expterm );
						F_B 				= func_output.F_value;
						f_eval++;

						if ( F_B - F_old <= alpha_B * suf_descent_1 ) {
						
							x_linesearch = y_B;
							F = F_B;
							/* set x, and grad_f */
							for( i = 0; i < n; i++ ){
								step[i] = x_linesearch[i] - x[i];
								norm_step += step[i] * step[i];
								x[i] = x_linesearch[i];
								grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
							}
							grad_eval++;
							norm_step = sqrt( norm_step );
							break;

						} else {
							
							logistic_setExpTerm( x_linesearch, wTx, sigmoid, ysigmoid, expterm );
							/* set x, and grad_f */
							for( i = 0; i < n; i++ ){
								step[i] = x_linesearch[i] - x[i];
								norm_step += step[i] * step[i];
								x[i] = x_linesearch[i];
								grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
							}
							grad_eval++;
							norm_step = sqrt( norm_step );
							break;
						
						}

					}

				}


				if ( j > maxback ){
					/* set x, and grad_f */
					for( i = 0; i < n; i++ ){
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						
						x[i] = x_linesearch[i];
						grad_f[i] = -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
					
					}

					grad_eval++;
					norm_step = sqrt( norm_step );
					break;
				}

				/* while loop for line search ends */

				alpha 				*= xi;
				proj_output 		= project( x, d, alpha );
				nDifforthant 		= proj_output.nV;
				x_linesearch 		= proj_output.project_vector;
				sameOrthant_prev 	= sameOrthant;
				sameOrthant 		= proj_output.same_sign;
				logistic_setExpTerm( x_linesearch, wTx, sigmoid, ysigmoid, expterm );
				func_output 		= logistic_func( x_linesearch, expterm );
				F 					= func_output.F_value;
				f 					= func_output.f_value;
				x_l1_norm 			= func_output.x_l1_norm;
				f_eval++;
				j++;


			}	
			

		} else {


			beta_iter++;

			select_beta();

			iter_type 			= 0;
			TRbeta_norm_beta 	= -( TRradiusBeta / norm_Beta );

			sprintf( type_iteration, "%s", "beta" );
			sub_iters 			= 1;


			for ( i = 0; i < n; i++ ){
				d[i] = TRbeta_norm_beta * Beta[i];
				x_linesearch[i] = x[i] + d[i];
				if( d[i] > 0 ){
					grad_F[i] = grad_f[i] + lambda;
				}else if( d[i] < 0 ){
					grad_F[i] = grad_f[i] - lambda;
				}else{
					grad_F[i] = 0.0;
				}
			}

			norm_d 				= TRradiusBeta;
			/* Perform line search to get updated x. */
			j 					= 0;
			alpha 				= 1.0;
			dirDer 				= dot_n5( grad_F, d, n );
			F_old 				= F;
			logistic_setExpTerm( x_linesearch, wTx, sigmoid, ysigmoid, expterm );
			func_output 		= logistic_func( x_linesearch, expterm );
			F 					= func_output.F_value;
			f 					= func_output.f_value;
			x_l1_norm 			= func_output.x_l1_norm;	
			f_eval++;

			suf_descent_1 		= param.eta * dirDer; 
			tmp 				= 0.0; /* used for update grad_f */
			
			// fprintf( stdout, "line search start\n" );

			/* perform line search */
			while( 1 ){
				
				// fprintf( stdout, "j = %d, dirDer: %f, F - F_old: %f\n", j, dirDer, F - F_old );
				if( F - F_old <= alpha * suf_descent_1 ){
					/* set x, and grad_f */
					for( i = 0; i < n; i++ ){
						step[i] 	= x_linesearch[i] - x[i];
						norm_step 	+= step[i] * step[i];
						x[i] 		= x_linesearch[i];
						grad_f[i] 	= -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
					}
					grad_eval++;
					norm_step 		= sqrt( norm_step );
					sprintf( linesearch_flag, "%s", "descent" );
					break;
				}

				if ( j > maxback ){
					/* set x, and grad_f */
					for( i = 0; i < n; i++ ){
						step[i] 	= x_linesearch[i] - x[i];
						norm_step 	+= step[i] * step[i];
						x[i] 		= x_linesearch[i];
						grad_f[i] 	= -1.0 / dm * dot_r( ysigmoid, problem.cX[i] );
					}
					grad_eval++;
					norm_step 		= sqrt( norm_step );
					break;
				}
				j++;
				alpha 		*= xi;	

				/* set new trial step */
				for( i = 0 ; i < n; i++ ){
					x_linesearch[i] = x[i] + alpha * d[i];
				}

				logistic_setExpTerm( x_linesearch, wTx, sigmoid, ysigmoid, expterm );
				func_output = logistic_func( x_linesearch, expterm );
				F 			= func_output.F_value;
				f 			= func_output.f_value;
				x_l1_norm 	= func_output.x_l1_norm;
				f_eval++;
			}		

			// printe( x, n );
			

		} /* end beta step */

		/* Keep running total of sub-problem iterations */

			
		/* output based on different iteration type */
		if ( param.print_level >= 4 ) {

			fo = fopen( param.output_file, "ab" );

			if ( !strcmp( type_iteration, "phi" ) ) {
				// fprintf(stdout, "outputting \n" );
				fprintf( fo, " %4s %5d %7.1e %7s %6d %5d %8.2e %8.2e %5d %5d %8.2e |", type_iteration, nI_phi, TRradiusPhi, output_cp.sub_prob_flag, rsType, sub_iters, sp_residual, sp_target_residual, output_cp.nV, nVmax, norm_d );
				fprintf( fo, " %3d %8.2e %7s %5d |\n", j, alpha, linesearch_flag, nDifforthant );
				// %4d %4d %1.3f %1.3f | %1.3f %1.3f %1.3f\n",  j, alpha,  );
			} else if ( !strcmp( type_iteration, "beta" ) ) {
				fprintf( fo, " %4s %5d %7.1e ------- %6d %5d -------- -------- ----- ----- %8.2e |", type_iteration, nI_phi, TRradiusBeta, rsType, sub_iters, norm_d );
				fprintf( fo, " %3d %8.2e %7s %5d |\n", j, alpha, linesearch_flag, 0 );
			}

			fclose( fo );
		}


		for ( i = 0; i < n; i++ ) {
			d[i] = 0.0;
		}

		nI_phi 	= 0;
		nI_beta	= 0;
		
		/* Set trust-region radius for next iteration. */
		if( iter_type == 1 ){
			// TRradiusPhi 	= MAX( 1E-3, MIN( 1E3, 10 * norm_step ) );
			TRradiusPhi 	= MAX( 1E-1, MIN( 1E3, 10 * norm_step ) );
		}else{
			TRradiusBeta 	= MAX( 1E-5, MIN( 1.0, norm_step ) );
		}
		
		norm_step = 0.0;

		setBetaPhi();
		iter++;


	}

}

Func_output_t logistic_func( double *x, double *expterm ){
	
	Func_output_t 	func_output;
 	int 			i;
   
	func_output.f_value = 0.0;
	for ( i = 0; i < m; i++ ) {
		func_output.f_value += log( 1.0 + expterm[i] );
	}

	func_output.f_value 	/= (double) m;
	func_output.x_l1_norm 	= l1_n5( x, n );
	func_output.F_value 	= func_output.f_value + lambda * func_output.x_l1_norm;

	return func_output;

}

void logistic_setExpTerm( double *x, double *wTx, double *sigmoid, double *ysigmoid, double *expterm ){
	
	int i;
   	
   	for( i = 0; i < m; i++ ) {
    
    	wTx[i] 		 = dot_c( x, problem.rX[i] );
		expterm[i] 	 = exp( -1.0 * y[i] * wTx[i] );
      	sigmoid[i] 	 = expterm[i] / ( 1.0 + expterm[i] );
      	ysigmoid[i]  = y[i] * sigmoid[i]; 
   		diagonal[i]  = sigmoid[i] * ( 1.0 - sigmoid[i] );
   }

}


double *logistic_hessVecProd( double *v ){

	double 	*hv;
	double 	*tmp_vector;
	int 	i;

	hv 			= Malloc( double, nI_phi_selected );
	tmp_vector 	= Malloc( double, m );


	for( i = 0; i < m; i++ ) {
		tmp_vector[i] = diagonal[i] * dot_c( v, problem.rXS[i] ); 
	}

	for( i = 0; i < nI_phi_selected; i++ ) {
		hv[i] = 1.0 / (double) m * dot_r( tmp_vector, problem.cXS[i] ) + 1E-8 * v[i];				
	}

	free (tmp_vector);
	return hv;

}


void logistic_setXF(){
	
	int 	i;
	int 	nnz;
	
	
	nnz = 0;

	/* set column format data */
	problem.cXS = Malloc( Node_rc_t *, nI_phi_selected );

	for(  i = 0; i < nI_phi_selected; i++ ){
		problem.cXS[i] = &X_col[col_ptr_head[I_phi_selected[i]]];
	}

	problem.rXS = Malloc( Node_rc_t *, m );
	
    for( i = 0; i < m + 1; i++ ){
        row_ptr_sp[i] = 0;
    }

    for( i = 0; i < nI_phi_selected; i++ ){
        Node_rc_t *x = problem.cXS[i];
        while( x->r_index != -1 ){
            nnz++;
            row_ptr_sp[x->r_index+1 ]++;
            x++;
        }
    }
    
    for( i = 1; i < m + 1; i++ ){
        row_ptr_sp[i] += row_ptr_sp[i-1] + 1;
    }

    X_row_S = Malloc( Node_rc_t, nnz + m );

    for( i = 0; i < m; i++ ){
        problem.rXS[i] = &X_row_S[row_ptr_sp[i]];
    }

    for( i = 0; i < nI_phi_selected; i++ ){
        Node_rc_t *x = problem.cXS[i];
        while( x->r_index != -1 ){
            int idx = x->r_index; 
            X_row_S[row_ptr_sp[idx]].c_index = i; 
            X_row_S[row_ptr_sp[idx]].r_index = idx;
            X_row_S[row_ptr_sp[idx]].value = x->value;
            row_ptr_sp[idx]++;
            x++;
        }
    }

    for( i = 0; i < m; i++ ){
        X_row_S[row_ptr_sp[i]].c_index = -1;
        X_row_S[row_ptr_sp[i]].r_index = -1;
    }

}

/* **********  Logistic Regression Loss Section End ******** */


/* **********  Least Square Loss Section Start ******** */



/* **********  Least Square Loss Section End ******** */

/* calculate sparsity in x */
double calculate_sparsity( double *x, int n ){

	int 	i;
	int 	zero_count;
	double 	zero_perc;

	zero_count = 0;

	for( i = 0; i < n; i++ ){
		if( x[i] == 0 ){
			zero_count++;
		}
	}

	zero_perc = (double) zero_count / (double) n;

	return zero_perc;
	
}


/* Sparse Operators Start */

double dot_r( const double *v, const Node_rc_t *x ){
   
   double 	vTx;
   
   vTx = 0.0;
   
   while( x->r_index != -1 ){
      vTx += v[x->r_index] * x->value;
      x++;
   }
   
   return vTx;

}

double dot_r_divide_num( const double *v, const Node_rc_t *x, const double num  ){
   
   double 	vTx;
   
   vTx = 0.0;
   
   while( x->r_index != -1 ){
      vTx += v[x->r_index] / num * x->value;
      x++;
   }
   
   return vTx;

}

double dot_c( const double *v, const Node_rc_t *x ){

   double vTx;

   vTx = 0.0;

   while( x->c_index != -1 ){
      vTx += v[x->c_index] * x->value;
      x++;
   }
   return vTx;
}


double norm2_sq_r( const Node_rc_t *x ){
   
   double val;

   val = 0.0;
   
   while( x->r_index != -1 ){
      val =+ x->value * x->value;
      x++;
   }
   return val;
}

double norm2_sq_c( const Node_rc_t *x ){
   
   double val;
   
   val = 0.0;

   while( x->c_index != -1 ){
      val =+ x->value * x->value;
      x++;
   }
   return val;
}

double dot_S_r( const double *v, const Node_rc_t *x, int *S ){
   
   double vTx;

   vTx = 0.0;

   while( x->r_index != -1 && *S != -1 ){
      if( x->r_index == *S ){
         vTx += *v * x->value;
         v++;
         x++;
         S++;
      }else if( x->r_index > *S ){
         S++;
         v++;
      }else if( x->r_index < *S ){
         x++;
      }
   }
   
   return vTx;  
}

double dot_S_c( const double *v, const Node_rc_t *x, int *S ){
   
   double vTx;
   vTx = 0.0;

   while( x->c_index != -1 && *S != -1 ){
      if( x->c_index == *S ){
         vTx += *v * x->value;
         v++;
         x++;
         S++;
      }else if( x->c_index > *S ){
         S++;
         v++;
      }else if( x->c_index < *S ){
         x++;
      }
   }
   
   return vTx;  
}

/* Sparse Operators End */

/* n5 operators start */
double dot_n5( double *v1, double *v2, int n ){
   
   int 		i, n5;
   double 	result;
   
   result 	= 0.0;

   if( n <= 0 ) {
   		return result;
   }

   n5 		= n % 5;

   for( i = 0; i < n5; i++ ){
      result += v1[i] * v2[i];
   }

   for( ; i < n; i += 5 ){
      result += v1[i]*v2[i] + v1[i+1]*v2[i+1] + v1[i+2]*v2[i+2] + v1[i+3]*v2[i+3] + v1[i+4]*v2[i+4];
   }

   return result;
}

double l1_n5( double *v, int n ){
   int 		i, n5;
   double 	result;

   result 	= 0.0;

   if( n <= 0 ) {
   		return result;
   } 

   n5 		= n % 5;

   for( i = 0; i < n5; i++ ){
      result += fabs(v[i]);
   }

   for( ; i < n; i += 5 ){
      result += fabs(v[i]) + fabs(v[i+1]) + fabs(v[i+2]) + fabs(v[i+3]) + fabs(v[i+4]);
   }
   return result;
}

/* n5 operators end */

/* print tools */
void print( double *vector, int n ){
   
   int i;
   for( i = 0; i < n; i++ ){
      fprintf(stdout, "%f ", vector[i] );
   }
   fprintf(stdout, "\n" );
}

void printInt( int *vector, int n ){
   int i;
   for( i = 0; i < n; i++ ){
      fprintf(stdout, "%d ", vector[i] );
   }
   fprintf(stdout, "\n" );
} 

void printe( double *vector, int n ) {

   int i;
   for( i = 0; i < n; i++ ){
      fprintf(stdout, "%e ", vector[i] );
   }
   fprintf(stdout, "\n" );

}

/* reduced space sub-problem solver */
Output_sp_t Directsolver( Input_direct_t input_direct ){

	Output_sp_t output_direct;



	return output_direct;

}


/* CG reduced-space solver */
Output_sp_t CGsolver( Input_CG_t input_cg ){

	Output_sp_t output_cg;
	
	int 		i;
	int 		nV;
	int 		iter_cg;
	int 		*I;
	int 		nI;
	int 		n;
	int 		nVmaxAllow;
	int 		max_CG_iter;
	int 		max_loop;
	double 		*x;
	double 		*grad_F;
	double 		TRradius;
	double 		*d_full;
	double 		*d_reduced;
	double 		*p;
	double 		*r;
	double 		*Hp;
	double 		norm_r0;
	double 		norm_r;
	double 		normrSq;
	double 		res_target;
	double 		pTHp;
	double 		alphaCG;
	double 		betaCG;
	double 		norm_d;
	double 		norm_old;
	char		sub_prob_flag[10];

	// for test
	double 		norm_x;

	nV 			= 0;
	iter_cg 	= 0;
	I 			= input_cg.I;
	nI 			= input_cg.nI;
	n 			= input_cg.n;
	nVmaxAllow 	= input_cg.nVmaxAllow;
	max_CG_iter = input_cg.maxCG_iter;
	max_loop 	= MIN( max_CG_iter,nI );
	x 			= input_cg.x;
	grad_F 		= input_cg.grad_F;
	TRradius 	= input_cg.TRradius;
	d_full 		= Malloc( double, n );
	d_reduced 	= Malloc( double, nI );
	p 			= Malloc( double, nI );
	r 			= Malloc( double, nI );
	Hp 			= Malloc( double, nI );
	norm_r0 	= 0.0;
	res_target 	= 0.0;
	pTHp 		= 0.0;
	alphaCG 	= 0.0;
	betaCG 		= 0.0;
	norm_d 		= 0.0;
	norm_old 	= 0.0;


	for( i = 0; i < nI; i++ ){
		r[i] 	= grad_F[ I[i] ];
		p[i] 	= -grad_F[ I[i] ];
		norm_r0 += p[i] * p[i];
		d_reduced[i] = 0.0;
	}
	

	for( i = 0; i < n; i++ ){
		d_full[i] = 0.0;
		norm_x += x[i] * x[i];
	}

	norm_x = sqrt( norm_x );

	normrSq = norm_r0;
	norm_r0 = sqrt( norm_r0 );
	norm_r 	= norm_r0;

	if( input_cg.rsType == 1 ){
		res_target = MAX( input_cg.eta_r * norm_r0, 1E-12 );
	}else{
		res_target = MAX( MIN( input_cg.eta_r, norm_r0 ) * norm_r0, 1E-12 );
	}
	
	while(1){

		/* Compute next linear CG iterate. */
		Hp = input_cg.hessVecProd(p);
		hessVec_eval++;
		pTHp = dot_n5( p, Hp, nI );
		alphaCG = normrSq / pTHp;

		norm_old = norm_r;
		norm_r = 0.0;
		nV = 0;
		norm_d = 0.0;
		iter_cg++;

		for( i = 0; i < nI; i++ ){
			d_reduced[i] += alphaCG * p[i];
			norm_d += d_reduced[i] * d_reduced[i];
			d_full[I[i]] = d_reduced[i];
			r[i] = r[i] + alphaCG * Hp[i];
			norm_r += r[i] * r[i];
		}

		normrSq = norm_r;
		norm_d = sqrt( norm_d );
		norm_r = sqrt( norm_r );

		/* calculate violating sets */
		for( i = 0; i < n; i++ ){
			if( x[i] > 0 && x[i] + d_full[i] < 0 ){
				nV++;
			}else if( x[i] < 0 && x[i] + d_full[i] > 0){
				nV++;
			}
		}

		/* Check for termination of CG. */
		if( norm_r <= res_target ){
			strncpy( sub_prob_flag, "tol:", sizeof(sub_prob_flag)-1);
			/* fprintf(stdout, "CG teriminate : CG residual\n" ); */
			break;
		}else if( nV > nVmaxAllow ){
			strncpy( sub_prob_flag, "vio:", sizeof(sub_prob_flag)-1);
			/* fprintf(stdout, "CG teriminate : CG violate\n" ); */
			break;
		}else if( norm_d >= TRradius ){
			strncpy( sub_prob_flag, "big:", sizeof(sub_prob_flag)-1);
			/* fprintf(stdout, "CG teriminate : CG big\n" ); */
			break;
		}else if( iter_cg > max_loop ){
			strncpy( sub_prob_flag, "max:", sizeof(sub_prob_flag)-1);
			/* fprintf(stdout, "CG teriminate : CG maximum\n" ); */
			break;
		}

		betaCG = normrSq / ( norm_old * norm_old );

		for( i = 0; i < nI; i++ ){
			p[i] = -r[i] + betaCG * p[i];
		}

	}

	/* calculate dirder */
	output_cg.dir_der = 0.0;
	for( i = 0; i < nI; i++ ){
		output_cg.dir_der += d_reduced[i] * grad_F[I[i]];
	}
	output_cg.d_full 		= d_full;
	output_cg.nV 			= nV;
	output_cg.iter 			= iter_cg;
	output_cg.res 			= norm_r;
	output_cg.norm_d 		= norm_d;
	output_cg.res_target 	= res_target;
	output_cg.nVmax 		= nVmaxAllow;
	sprintf( output_cg.sub_prob_flag, "%s", sub_prob_flag );

	return output_cg;	

}

/* CD reduced-space solver */
/*
	A combination of a standard coordiante descent (CD) method, followed by
	a search along the dogleg path that is built upon the Cauchy point and
	vector resulting from the CD minimization procedure.  Only works on a 
	subset of the data that only takes takes columns corresponding to I.
*/
Output_sp_t CDsolver_DL_diffd_logloss( Input_CD_t input_cd ){

	// fprintf( stdout, "Use CD solver\n" );
	
	Output_sp_t output_cd;

	int 		*I;
	int 		*V_p;
	int 		*V_n;
	int 		*V;
	int 		*x_I_pos;
	int 		*x_I_neg;
	int 		nV;
	int 		nI;
	int 		n;
	int 		m;
	int 		i;
	int 		j;
	int 		max_CD_iter;
	int 		iter_cd;
	int 		nVmax;
	int 		HvProds;
	int 		back_tracks;
	int 		max_back_tracks;
	int 		cur_row;
	int 		terminate;
	double 		*grad_F;
	double 		*x_I;
	double 		*D;
	double 		*H_diags; /* diagonal of hessian matrix 1/m * X^TDX */
	double 		*g;
	double 		*x;
	double 		*d_r; /* reference direction */
	double 		*d_full;
	double 		*d_I;
	double 		*d_I_old;
	double 		*X_I_d;
	double 		*p;
	double 		*Hg;
	double 		*temp;
	double 		*diff;
	double 		a;
	double 		b;
	double 		c;
	double 		G_j;
	double 		tmp;
	double 		mu;
	double 		norm_g;
	double 		norm_g_sq;
	double 		norm_d_r;
	double 		alpha;
	double 		beta;
	double 		norm_r_sq;
	double 		norm_r;
	double 		norm_d;
	double 		norm_d_I;
	double 		eta_r;
	double 		res;
	double 		res_target;
	double 		frac_viol;
	double 		gTHg;
	double 		descent;
	double 		suf_descent;
	double 		TRradius;
	double 		nVmaxAllow;
	char		sub_prob_flag[10];

	Node_rc_t 	*X_col_j; /* data matrix column */

	I 			= input_cd.I;
	nI 			= input_cd.nI;
	n 			= input_cd.n;
	m 			= input_cd.m;
	max_CD_iter	= input_cd.max_CD_iter;
	nVmaxAllow 	= input_cd.nVmaxAllow;
	frac_viol 	= input_cd.fracViol;
	TRradius 	= input_cd.TRradius;
	D 			= input_cd.help_vector;
	x_I 		= input_cd.x_I;
	x_I_pos 	= input_cd.x_I_pos;
	x_I_neg 	= input_cd.x_I_neg;
	eta_r 		= input_cd.eta_r;
	grad_F 		= input_cd.grad_F;
	V 			= Malloc( int, nI );
	V_p 		= Malloc( int, nI );
	V_n 		= Malloc( int, nI );
	d_full		= Malloc( double, n );
	d_r 		= Malloc( double, nI );
	d_I 		= Malloc( double, nI );
	d_I_old 	= Malloc( double, nI );
	p 			= Malloc( double, nI );
	g 			= Malloc( double, nI );
	Hg 			= Malloc( double, nI );
	H_diags		= Malloc( double, n );
	X_I_d 		= Malloc( double, m );
	temp 		= Malloc( double, m );

	res_target 	= 0.0;
	res 		= 0.0;
	norm_d_I 	= 0.0;
	norm_d 		= 0.0;
	terminate 	= FALSE;
	HvProds 	= 0;
	norm_g_sq 	= 0.0;

	/* initialize declared vectors */
	for ( i = 0; i < n; i++ ) {
		H_diags[i] 	= 0.0;
		d_full[i]	= 0.0;
	}

	for ( i = 0; i < m; i++ ) {
		X_I_d[i] 	= 0.0;
		temp[i] 	= 0.0;
	}

	for ( i = 0; i < nI; i++ ) {

		V[i] 		= 0;
		V_p[i] 		= 0;
		V_n[i]		= 0;
		d_r[i] 		= 0.0;
		d_I[i] 		= 0.0;
		d_I_old[i] 	= 0.0;
		p[i] 		= 0.0;
		g[i] 		= 0.0;
		Hg[i] 		= 0.0;

	}


	// fprintf(stdout, "x_I: " );
	// print( x_I, nI );
	// fprintf(stderr, "grad_F:" );
	// print( grad_F, n );
	

	for( i = 0; i < nI; i++ ){

		g[i] 	= grad_F[ I[i] ];
		p[i] 	= -grad_F[ I[i] ];
		norm_g_sq += g[i] * g[i];
		d_I[i] 	= 0.0;
	
	}

	norm_g 		= sqrt( norm_g_sq );
	nVmax 		= ceil( MIN( nVmaxAllow, frac_viol * nI ) );
	// fprintf(stdout, "nVmaxAllow: %f, frac_viol: %f, nI: %d\n", nVmaxAllow, frac_viol, nI );
	Hg 			= input_cd.hessVecProd( g );
	hessVec_eval++;
	HvProds 	= 1;
	nV 			= 0;
	gTHg 		= dot_n5( g, Hg, nI );

	/* Minimizer of model along -alpha*g. Cauchy point. */
	alpha 		= norm_g_sq / gTHg;

	for( i = 0; i < nI; i++ ){
		d_r[i] 	= -alpha * g[i];
	}
	norm_d_r 	= alpha * norm_g;

	suf_descent = - alpha * norm_g_sq;


	for( i = 0; i < n; i++ ){
		d_full[i] = 0.0;
	}

	/* Set weak (1) or tight (2) termination depending on value of rsType. */
	if( input_cd.rsType == 1 ){
		res_target = MAX( input_cd.eta_r * norm_g, 1E-12 );
	} else {
		res_target = MAX( MIN( input_cd.eta_r, norm_g ) * norm_g, 1E-12 );
	} 

	// fprintf(stdout, "res_target: %f\n", res_target );
	// getchar();
	/* 
		Reference direction is larger than trust-region radius. Perform 
    	backtracking along the reference direction in order to satisfy that
    	the number of variables (nV) that switch orthants is less than the
    	maximum allowed quantity (nVmax) computed above, and that the step
    	satisfies the implicit trust-region constraint.
    */
	if ( TRradius <= norm_d_r ) {

		// fprintf(stdout, "case 1, TRradius: %f, norm_d_r: %f\n", TRradius, norm_d_r );

		back_tracks = 0;
		beta 		= ( TRradius / norm_g );

 		// fprintf(stderr, "beta: %f\n", beta );

 		while( terminate == FALSE ) {

			for ( i = 0; i < nI; i++ ) {
				d_I[i] 		= - beta * g[i];
				if ( x_I[i] + d_I[i] < 0.0 && x_I_pos[i] == 1 ) {
					V_p[i] 	= 1;
					V[i] 	= 1;
					nV++;
				} else if ( x_I[i] + d_I[i] > 0.0 && x_I_neg[i] == 1 ) {
					V_n[i] 	= 1;
					V[i] 	= 1;
					nV++;
				}
			}

			if ( nV <= nVmax ) {
				norm_d_I 	= beta * norm_g;
				sprintf( sub_prob_flag, ":dR%d", back_tracks );
				terminate 	= TRUE;
			} else {
				back_tracks++;
				beta *= 0.5;
			}

			for ( i = 0; i < nI; i++ ) {
				V_p[i] 	= 0;
				V_n[i] 	= 0;
				V[i] 	= 0;
			}
			nV 		= 0;

 		}

 		res 		= 0.0;
 		for ( i = 0; i < nI; i++ ) {
 			res 	+= ( g[i] - beta * Hg[i] ) * ( g[i] - beta * Hg[i] );
 		}
 		res 		= sqrt( res );
 		iter_cd 	= 0;
 		
	} 
	/*
    	The reference direction dR is smaller than the trust-region radius.
    	Perform coordinate descent for an approximate solution and use it.
	*/
	else {
		// fprintf(stdout, "case 2, TRradius: %f, norm_d_r: %f\n", TRradius, norm_d_r );
		copy( d_I_old, d_I, nI );
		iter_cd 	= 1;

		/* Main CD loop */
		while ( terminate == FALSE ) {

			/* Shuffle first nI integers for random coordinate descent */
			// for ( i = 0; i < nI; i++ ) {
			// 	j 	= i + rand() % ( nI - i );
			// 	swap( &I[i], &I[j] );
			// }

			for ( i = 0; i < nI; i++ ) {

				j 			= I[i];
				H_diags[j] 	= mu;	

				X_col_j 	= problem.cXS[i];
				tmp 		= 0.0;

				/* calculate j-th diagonal entry in Hessian matrix */
				/* need to check whether do single time or multiple times */
				// double norm_temp = 0.0; // for test
				while( X_col_j->r_index != -1 ){
					H_diags[j] += X_col_j->value * X_col_j->value * D[X_col_j->r_index] / (double) m;
					temp[ X_col_j->r_index ] = D[ X_col_j->r_index ] * X_col_j->value;
					// norm_temp += temp[ X_col_j->r_index ] * temp[ X_col_j->r_index ];
					X_col_j++;
					
				}
				// fprintf(stderr, "line 1742, g[i]: %f\n", g[i] );
				// print( X_I_d, m );
				G_j 		= g[i] + dot_n5( temp, X_I_d, m ) / (double) m ;

				alpha 		= - G_j / H_diags[j];

				/* update the j-th entry */
				d_I[i] 		+= alpha;
				/* update X_I_d */
				X_col_j 	= problem.cXS[i];
				
				while( X_col_j->r_index != -1 ) {
					X_I_d[ X_col_j->r_index ] += alpha * X_col_j->value;
					X_col_j++;
				}
				// fprintf(stdout, "norm temp: %f, g[i]: %f\n", sqrt( norm_temp ), g[i] );
				// fprintf(stdout, "j: %d, Hjj: %f, G_j: %f, alpha: %f, d_I[i]: %f\n", j, H_diags[j], G_j, alpha, d_I[i] );
				
				/* reset temp */
				for ( j = 0; j < m; j++ ) {
					temp[j] 	= 0.0;
				}
			}

			/* Norm of step computed so far */
			norm_r 		= 0.0;
			norm_d_I 	= 0.0;

			for ( i = 0; i < nI; i++ ) {	
				norm_d_I 	+= d_I[i] * d_I[i];
				norm_r 		+= ( d_I[i] - d_I_old[i] ) * ( d_I[i] - d_I_old[i] );
				if ( x_I_pos[i] == 1 && x_I[i] + d_I[i] < 0.0 ) {
					V_p[i] 	= 1;
					V[i] 	= 1;
					nV++;
				} else if ( x_I_neg[i] == 1 && x_I[i] + d_I[i] > 0.0 ) {
					V_n[i] 	= 1;
					V[i] 	= 1;
					nV++;
				}
			}

			norm_d_I 		= sqrt( norm_d_I );
			// fprintf(stderr, "norm_d_I: %f\n", norm_d_I );
			// getchar();
			norm_r 			= sqrt( norm_r );
			/* Descent associated with current estimate dI. */
			descent  		= dot_n5( g, d_I, nI );

			// fprintf(stdout, "nV: %d, descent: %f, nVmax: %d, suf_descent: %f\n", nV, descent, nVmax, suf_descent );
			// printInt( V_p, nI );
			// printInt( V_n, nI );
			// printInt( V, nI );

			if ( descent <= suf_descent && norm_r <= res_target ) {
				strncpy( sub_prob_flag, "tol:", sizeof(sub_prob_flag)-1);
				terminate 	= TRUE;

			} else if ( descent <= suf_descent && nV > nVmax ) {
				strncpy( sub_prob_flag, "vio:", sizeof(sub_prob_flag)-1);
				terminate 	= TRUE;
			} else if ( descent <= suf_descent && norm_d_I >= TRradius ) {
				strncpy( sub_prob_flag, "big:", sizeof(sub_prob_flag)-1);
				terminate 	= TRUE;
			} else if ( iter_cd >= max_CD_iter ) {
				strncpy( sub_prob_flag, "max:", sizeof(sub_prob_flag)-1);
				terminate 	= TRUE;
			} else {
				terminate 	= FALSE;
			}

			if ( terminate == TRUE ) {
				// fprintf(stdout, "termination flag: %s\n", sub_prob_flag );
				// fprintf(stdout, "max_CD_iter: %d\n", max_CD_iter );
				res 		= norm_r;
				break;
			}

			/* Save current solution estimate. */
			copy( d_I_old, d_I, nI );
			
			/* Update the counter. */
			iter_cd++;

			/* reset some vector */
			for ( i = 0; i < nI; i++ ) {
				V_p[i] 	= 0.0;
				V_n[i] 	= 0.0;
				V[i]	= 0.0;
			}

			for ( i = 0; i < m; i++ ) {
				temp[i] = 0.0;
			}

			for ( i = 0; i < n; i++ ) {
				H_diags[i] = 0.0;
			}
			nV 	 = 0;
	
		} /* end of while(TRUE) Main CD loop */

		diff 	= Malloc( double, nI );
		a 		= 0.0;
		b 		= 0.0;
		c 		= norm_d_r * norm_d_r - TRradius * TRradius;

 	   /* 
 	   		Use the reference direction dR and the just computed CD direction to 
 		   	form the dogleg path.  Do backtracking along this path to ensure
    		that the descent condition and the number of variables (nV) that
    		switch orthants is less than the quantity nVmax computed above.
		*/

		for ( i = 0; i < nI; i++ ) {
			diff[i] 	= d_I[i] - d_r[i];
			a 			+= diff[i] * diff[i];
			b 			+= 2 * d_r[i] * diff[i];
		}

		beta 	= MIN( 1.0, -b + sqrt( b * b - 4 * a * c ) / ( 2 * a ) );

		back_tracks 	= 0;
		max_back_tracks = 5;
		descent 		= 0.0;
		nV 				= 0;
		for ( i = 0; i < nI; i++ ) {
			V_p[i] 	= 0.0;
			V_n[i] 	= 0.0;
			V[i]	= 0.0;
		}
		terminate 		= TRUE;
		
		// fprintf(stdout, "=========================\n" );
		// fprintf(stdout, "Search along dogleg path \n" );	
		// getchar();

		while( TRUE ) {

			norm_d_I 		= 0.0;

			for ( i = 0; i < nI; i++ ) {
				d_I[i] 		= d_r[i] + beta * diff[i];
				norm_d_I 	+= d_I[i] * d_I[i];
				descent 	+= g[i] * d_I[i];
				if ( x_I[i] + d_I[i] < 0.0 && x_I_pos[i] == 1 ) {
					V_p[i] 	= 1;
					V[i] 	= 1;
					nV++;
				} else if ( x_I[i] + d_I[i] > 0.0 && x_I_neg[i] == 1 ) {
					V_n[i] 	= 1;
					V[i] 	= 1;
					nV++;
				}
			}
			// fprintf(stdout, "nV: %d, descent: %f, nVmax: %d, suf_descent: %f\n", nV, descent, nVmax, suf_descent );
			if ( nV > nVmax || descent > suf_descent ) {
				beta *= 0.5;
				back_tracks++;
			} else {
				norm_d_I 	= sqrt( norm_d_I );
				if ( beta == 1.0 ) {
					sprintf( sub_prob_flag, "%sCD", sub_prob_flag ); 
				} else {
					sprintf( sub_prob_flag, "%sDL%d", sub_prob_flag, back_tracks );
				}
				break;
			}

			// fprintf(stderr, "max_back_tracks: %d\n", max_back_tracks );
			if ( back_tracks > max_back_tracks ) {
				/*
					Searching on the segment [dR,dI] failed.  Now just
					perform backtracking along the reference direction dR.
				*/
				back_tracks = 0;
				beta 		= 1.0;
				V_p 		= Malloc( int, nI );
				V_n 		= Malloc( int, nI );
				V 			= Malloc( int, nI );
				nV 			= 0;

				while( TRUE ) {

					for ( i = 0; i < nI; i++ ) {
						d_I[i] = beta * d_r[i];
						if ( x_I[i] + d_I[i] < 0.0 && x_I_pos[i] == 1 ) {
							V_p[i] 	= 1;
							V[i] 	= 1;
							nV++;
						} else if ( x_I[i] + d_I[i] > 0.0 && x_I_neg[i] == 1 ) {
							V_n[i] 	= 1;
							V[i] 	= 1;
							nV++;
						}					
					}

					if ( nV <= nVmax ) {
						norm_d_I 	= beta * norm_d_r;
						sprintf( sub_prob_flag, "%sdR%d", sub_prob_flag, back_tracks ); 
						break;
					} else {
						back_tracks++;
						beta *= 0.5;
					}

					for ( i = 0; i < nI; i++ ) {
						V_p[i] 	= 0.0;
						V_n[i] 	= 0.0;
						V[i]	= 0.0;
					}
					nV 			= 0;

				}

				break;

			}

			/* reset */
			descent 	= 0.0;
			// V_p 		= Malloc( int, nI );
			// V_n 		= Malloc( int, nI );
			// V 			= Malloc( int, nI );
			nV 			= 0;			
			for ( i = 0; i < nI; i++ ) {
				V_p[i] 	= 0.0;
				V_n[i] 	= 0.0;
				V[i]	= 0.0;
			}
		}		
	
	} /* end of else */

 
	// fprintf(stdout, "nV: %d, descent: %f, nVmax: %d, suf_descent: %f\n", nV, descent, nVmax, suf_descent );


	output_cd.dir_der = 0.0;
	for ( i = 0; i < nI; i++ ) {
		output_cd.dir_der += d_I[i] * grad_F[ I[i] ];
		d_full[ I[i] ] = d_I[i];
	}


	output_cd.d_full 		= d_full;
	output_cd.nV 			= nV;
	output_cd.iter 			= iter_cd;
	output_cd.res 			= res;
	output_cd.norm_d 		= norm_d_I;
	output_cd.res_target 	= res_target;
	output_cd.nVmax 		= nVmax;
	sprintf( output_cd.sub_prob_flag, "%s", sub_prob_flag );
	
	// fprintf(stderr, "d_full: " );
	// print( d_full, n );
	// getchar();
	return output_cd;

}

/* CD reduced-space solver */
/*
	Compute the search direction based on a basic dog-leg search.  
*/
Output_CD_t CD_direct_solver( Input_CD_t input_cd ){
	
	Output_CD_t output_cd;

	int 		*S;
	int 		*I;
	int 		nI;
	int 		n;
	int 		i;
	int 		max_CD_iter;
	int 		iter_cd;
	int 		nV;
	double 		*x;
	double 		nVmaxAllow;
	double 		*r;
	double 		*d_full;
	double 		*d_reduced;
	double 		*p;
	double 		norm_r0;
	double 		norm_r;
	double 		norm_r_sq;
	double 		norm_d;
	double 		eta_r;
	double 		res_target;

	I 			= input_cd.I;
	nI 			= input_cd.nI;
	n 			= input_cd.n;
	nVmaxAllow 	= input_cd.nVmaxAllow;
	d_full		= Malloc( double, n );
	d_reduced 	= Malloc( double, nI );
	p 			= Malloc( double, nI );
	r 			= Malloc( double, nI );

	norm_r0 	= 0.0;
	res_target 	= 0.0;
	norm_d 		= 0.0;

	for( i = 0; i < nI; i++ ){
		r[i] = grad_F[ S[i] ];
		p[i] = -grad_F[ S[i] ];
		norm_r0 += p[i] * p[i];
		d_reduced[i] = 0.0;
	}


	for( i = 0; i < n; i++ ){
		d_full[i] = 0.0;
	}

	if( input_cd.rsType == 1 ){
		res_target = MAX( input_cd.eta_r * norm_r0, 1E-10 );
	} else {
		res_target = MAX( MIN( input_cd.eta_r, norm_r0 ) * norm_r0, 1E-10 );
	} 


	norm_r_sq  	= norm_r0;
	norm_r0 	= sqrt( norm_r0 );
	norm_r 		= norm_r0;

	while( TRUE ){

		for( i = 0; i < nI; i++ ){
			


		}

		if( norm_r <= res_target ){

			break;


		} else if ( iter_cd > max_CD_iter ){

			break;
		
		}

		iter_cd++;
	}


	output_cd.dir_der = 0.0;

	for ( i = 0; i < nI; i++ ) {
		output_cd.dir_der += d_reduced[i] * grad_F[ S[i] ];
	}

	output_cd.d_full 	= d_full;
	output_cd.nV 		= nV;
	output_cd.iter 		= iter_cd;


	return output_cd;

}

/************** Read parses and settings ************/
void parse_command_line( int argc, char **argv ){

	int i;
	param.max_iter 		= 1000;
	param.print_level 	= 2;
	param.Gamma 		= 1.0;
	param.eta_r 		= 1E-1;
	param.eta 			= 1E-2;
	param.xi 			= 0.5;
	param.tol_absolute 	= 1E-6;
	param.tol_relative 	= 1E-6;
	param.betaFrac 		= 1.0;
	param.phiFrac 		= 1.0;
	param.tryCG 		= 0;
	param.tryCD 		= 1;
	param.tryCrossOver 	= 0;
	param.max_CG_iter 	= MAX_VALUE;
	param.max_CD_iter 	= 1000;
	param.crossOverTol 	= 1E-1;
	param.TRradiusPhi 	= 1E3;
	param.TRradiusBeta 	= 1E-1;
	param.maxback 		= 100;
	param.fracViol 		= 0.1;
	param.nVmaxAllow 	= MAX_VALUE;
	param.checkResEvery = 1;
	param.termType 		= 2;
	param.lambdaCoeff 	= 1.0;
	param.scaling 		= 0;
	param.objective_function_type = 1;
	param.lambda 		= -1.0;
	param.sub_solver 	= 1;
	param.print_every 	= 10;
	param.max_time 		= 3600;


    for ( i = 1; i < argc; i++ ){
        if( argv[i][0] != '-' ) break;
        if( ++i >= argc ){
        	perror("Invalid arguments.");
            exit(0);
        }
        switch( argv[i-1][1] ){
            case 'p':
                strncpy( param.profile_file, argv[i], sizeof(param.profile_file)-1 );
                break;
            case 't':
            	strncpy( param.profile_file, "test/FaRSA_t.profile", sizeof(param.profile_file)-1 );
            	break;
            case 'n':
            	strncpy( param.name, argv[i], sizeof(param.name)-1);
            	break;
            case 's':
            	param.sub_solver = atoi( argv[i] );
            	break;
        }
    }
}

/* read setting from profile */
void load_setting_profile(){

	char *attribute;
	char *value;
	FILE *fp;

	fp = fopen( param.profile_file, "r" );

	if( fp == NULL ){
		fprintf( stdout, "Can not open profile file...\n" );
		perror("profile_file: fopen");
		exit(1);
	}

	max_line_size = 1024; 
	line = Malloc( char, max_line_size );

	while( get_line(fp) != NULL ){
		attribute = strtok( line, ":");
		
		value = strtok( NULL, " \t\n");

		/* If current value is null, then continue */
		if( value == NULL ){
			continue;
		}
		
		if( !strcmp( attribute, "maximum_iteration" ) ){
			param.max_iter = atoi(value);
			continue;
		}else if( !strcmp( attribute, "lambda_coefficient" ) ){
			param.lambdaCoeff = atof(value);
			continue;
		}else if( !strcmp( attribute, "Gamma" ) ){
			param.Gamma = atof(value);
			continue;
		}else if( !strcmp( attribute, "eta" ) ){
			param.eta = atof(value);
			continue;
		}else if( !strcmp( attribute, "eta_r" ) ){
			param.eta_r = atof(value);
			continue;
		}else if( !strcmp( attribute, "xi" ) ){
			param.xi = atof(value);
			continue;
		}else if( !strcmp( attribute, "xi" ) ){
			param.xi = atof(value);
			continue;
		}else if( !strcmp( attribute, "tolerance" ) ){
			param.tol_relative = atof(value);
			param.tol_absolute = atof(value);
			continue;
		}else if( !strcmp( attribute, "beta_fraction" ) ){
			param.betaFrac = atof(value);
			param.betaFrac = param.betaFrac > 1.0 ? 1.0 : param.betaFrac;
			continue;
		}else if( !strcmp( attribute, "phi_fraction" ) ){
			param.phiFrac = atof(value);
			param.phiFrac = param.phiFrac > 1.0 ? 1.0 : param.phiFrac;
			continue;
		}else if( !strcmp( attribute, "tryCG" ) ){
			param.tryCG = atoi(value);
			continue;
		}else if( !strcmp( attribute, "tryCD" ) ){
			param.tryCD = atoi(value);
			continue;
		}else if( !strcmp( attribute, "tryCrossOver" ) ){
			param.tryCrossOver = atoi(value);
			continue;			
		}else if( !strcmp( attribute, "max_CG_iter") ){
			/* atoi does not work for XEY */
			param.max_CG_iter = (int) atof(value);
			continue;
		}else if( !strcmp( attribute, "max_CD_iter") ){
			param.max_CD_iter = atoi(value);
			continue;
		}else if( !strcmp( attribute, "crossOverTol") ){
			param.crossOverTol = atof(value);
			continue;
		}else if( !strcmp( attribute, "TRradiusPhi") ){
			param.TRradiusPhi = atof(value);
			continue;
		}else if( !strcmp( attribute, "TRradiusBeta") ){
			param.TRradiusBeta = atof(value);
			continue;
		}else if( !strcmp( attribute, "maxback") ){
			param.maxback = atof(value);
			continue;
		}else if( !strcmp( attribute, "frac_viol") ){
			param.fracViol = atof(value);
			param.fracViol = param.fracViol > 1.0 ? 1.0 : param.fracViol;
			continue;
		}else if( !strcmp( attribute, "nV_max_allow") ){
			param.nVmaxAllow =  (int) atof(value);
			continue;
		}else if( !strcmp( attribute, "checkResEvery") ){
			param.checkResEvery = atoi(value);
			continue;
		}else if( !strcmp( attribute, "term_type") ){
			param.termType = atoi(value);
			continue;
		}else if( !strcmp( attribute, "objective_function_type") ){
			param.objective_function_type = atoi(value);
			continue;
		}else if( !strcmp( attribute, "data_file") ){
			int len;
			len = (int) strlen( (char *)value );
			if( len > 1024 ){
				fprintf(stdout, "Input path is too long\n" );
				exit(-1);
			}
			strncpy( param.prob_name, (char *)value, sizeof(param.prob_name)-1 );
			continue;
		}else if( !strcmp( attribute, "data_format") ){
			int len;
			len = (int) strlen( value );
			if( len > 1024 ){
				fprintf(stdout, "Data format is invalid.\n" );
			}
			strncpy( param.data_format, value, sizeof(param.data_format)-1);
			continue;
		}else if( !strcmp( attribute, "lambda") ){
			param.lambda = atof( value );
			continue;
		}else if( !strcmp( attribute, "output_file") ){
			int len;
			len = (int) strlen( value );
			if( len > 1024 ){
				fprintf(stdout, "Output file is invalid.\n" );
			}			
			strncpy( param.output_file, value, sizeof(param.output_file)-1);
			continue;
		} else if ( !strcmp( attribute, "print_level") ) {
			param.print_level = atoi( value );
			continue;
		} else if ( !strcmp( attribute, "sub_problem_solver" ) ) {
			param.sub_solver = atoi( value );
			continue;
		} else if ( !strcmp( attribute, "max_time" ) ) {
			param.max_time = atof( value );
			continue;
		}
	}

	fclose(fp);
}

/* initialize rest parameters */
void initialize_rest_param( struct Input_FaRSA input_farsa  ){

	/* initialize parameters */
	/* objective function is not personalized, then read problem */
	if( param.objective_function_type != 0 ){

		read_problem( &problem );

		m 			= problem.num_samples;
		n 			= problem.num_features;
		y 			= Malloc( double, m );
		dm 			= (double) m;
		lambda 		= param.lambda < 0 ? param.lambdaCoeff / (double) problem.num_samples : param.lambda;
		
	}else if( param.objective_function_type == 0 ){

		if( input_farsa.n == 0 ){
			fprintf(stdout, "Current n is 0. Please set n.\n" );
			exit(1);
		}
		lambda 		= param.lambda;
		n 			= input_farsa.n;
		m 			= 1;
		dm 			= (double) m;

	}else{
		
		fprintf(stdout, "Loss type is invalid\n" );
		exit(1);
	
	}

	iter 			= 0;
	beta_iter 		= 0;
	phi_iter 		= 0;
	sub_iters		= 0;

	f_eval 			= 0;
	grad_eval 		= 0;
	hessVec_eval 	= 0;

	rsType 			= 1;
	norm_beta 		= 0.0;
	norm_phi 		= 0.0;
	norm_phi0 		= 0.0;
	norm_beta0 		= 0.0;

	nI_phi 			= 0;
	nI_beta 		= 0;


	beta 			= Malloc( double, n );
	Beta 			= Malloc( double, n );
	phi 			= Malloc( double, n );

	grad_f 			= Malloc( double, n );
	grad_F 			= Malloc( double, n );

	I_phi 			= Malloc( int, n + 1 );
	I_beta 			= Malloc( int, n + 1 );
	I_z 			= Malloc( int, n );
	I_n 			= Malloc( int, n );
	I_p 			= Malloc( int, n );
	x 				= Malloc( double, n );
	Gamma 			= param.Gamma;
	tol_absolute 	= param.tol_absolute;
	tol_relative 	= param.tol_relative;

	d 				= Malloc( double, n );
	TRradiusBeta 	= param.TRradiusBeta;
	TRradiusPhi 	= param.TRradiusPhi;
	alpha 			= 1.0;
	x_linesearch 	= Malloc( double, n );
	dirDer 			= 0.0;
	F_old 			= 0.0;  /* objective function value for last step */
	F 				= 0.0; /* objective function value on current step */
	step 			= Malloc( double, n ); /* difference between two steps for trust region update */
	norm_step 		= 0.0;
	maxback 		= param.maxback;
	xi 				= param.xi;
	sameOrthant 	= TRUE;

	row_ptr_sp 		= Malloc( int, m + 1 );

}

/* function for selecting suitable read data */
void read_problem( Problem_t *prob ){
	if( param.scaling == 0 ){
		read_sparse_problem( prob );
	}else if( param.scaling == 1 ){
		read_scaling_problem( prob );
	}else if( param.scaling == 2 ){
		read_scaling_problem( prob );
	}
}

/* get line function */
char* get_line( FILE *input ){

	int len;

	if( fgets( line, max_line_size, input ) == NULL )
		return NULL;

	/* if current max_line_size doesn't find the \n, then enlarge 
	 the max_line_size */
	while( strrchr( line, '\n' ) == NULL ){
		max_line_size *= 2;
		// line = Realloc( char, max_line_size );
		line = (char *) realloc( line, max_line_size );
		len = (int) strlen(line);
		if( fgets( line+len, max_line_size-len, input ) == NULL )
			break;
	}
	return line;

}

void read_sparse_problem( Problem_t *prob ){

	int max_index;
	char *label, *idx, *val;
	size_t entries;

	FILE *fp = fopen( param.prob_name, "r" );

    if( fp == NULL ){
        fprintf( stdout, "Can't open input file %s\n", param.prob_name );
        // exit(0) behave like return 0 in main() function, exit(1) behave like return 1 
        exit(1);
    }

	max_line_size = 1024;
	prob->num_samples = 0;
	entries = 0;

	
	line = Malloc( char, max_line_size );


	while( get_line(fp) != NULL ){
	
		label = strtok( line, " \t\n");
        // features
        while(1){

            idx = strtok( NULL, ":" );
            val = strtok( NULL, " \t" );

            if( val == NULL || idx == NULL )
                break;

            entries++;      
        }
        prob->num_samples++; 		
	}
	rewind( fp );



	prob->y = Malloc( double, prob->num_samples );
	prob->rX = Malloc( Node_rc_t *, prob->num_samples ); 

	/* sort by samples */
	X = Malloc( Node_rc_t, entries + prob->num_samples );

	max_index = -1;

	int j = 0;
	int i;

	for ( i = 0; i < prob->num_samples; i++ ){
	
		get_line(fp);

		prob->rX[i] = &X[j];

		/* strtok: The C library function char *strtok(char *str, const char *delim) breaks 
		 string str into a series of tokens using the delimiter delim. */
		label = strtok( line, " \t\n");
		
		if( label == NULL )
			continue;


		prob->y[i] = atof(label);

		while(1){
			/* The first call to strtok() returns a pointer to the first token in the string 
			 pointed to by s1. Subsequent calls to strtok() must pass a NULL pointer as the 
			 first argument, in order to get the next token in the string. */
			idx = strtok( NULL, ":" );
			val = strtok( NULL, " \t" );

			if( val == NULL || idx == NULL )
				break;

			X[j].c_index 	= atoi(idx) - 1; /* c_index is index of feature */
			X[j].r_index 	= i; /* r_index is index of sample */
			X[j].value 		= atof(val);
            
            if( max_index < X[j].c_index ){
                max_index = X[j].c_index;
            }
			j++;

		}	

		X[j].r_index 	= -1;
		X[j++].c_index 	= -1;

	}
	prob->num_features 	= max_index + 1;


	transpose( prob );

	if( param.print_level >= 1 && param.print_level != 3 && param.print_level != 4 ){
		fprintf(stdout, "Dataset Description:\n" );
		fprintf(stdout, "Number of samples : %d\n", problem.num_samples );
		fprintf(stdout, "Number of features: %d\n", problem.num_features );		
	}

	fclose(fp);

}


void transpose( struct Problem *prob ){

    int i;
    
    int local_m;
    int local_n;

    local_m = prob->num_samples;
    local_n = prob->num_features;
    
    int nnz;
    nnz = 0;


    col_ptr = Malloc( int, local_n + 1 );


	if ( col_ptr == NULL) {
		fprintf( stdout, "Memory issue, please wait for a few seconds and run again! \n" );
        exit(-1);
    }

    
    col_ptr_head 	= Malloc( int, local_n + 1 );

    prob->cX 		= Malloc( Node_rc_t *, prob->num_features );

    for( i = 0; i < local_n + 1; i++ ){
        col_ptr[i] = 0;
        col_ptr_head[i] = 0;
    }

    for( i = 0; i < local_m; i++ ){
        Node_rc_t *x = prob->rX[i];
        while( x->c_index != -1 ){
            nnz++;
            col_ptr[x->c_index+1 ]++;
            col_ptr_head[x->c_index+1]++;
            x++;
        }
    }

    
    for( i = 1; i < local_n + 1; i++ ){
        col_ptr[i] += col_ptr[i-1] + 1;
        col_ptr_head[i] += col_ptr_head[i-1] + 1;
    }


    X_col = Malloc( Node_rc_t, nnz + local_n );

    for( i = 0; i < local_n; i++ ){
        problem.cX[i] = &X_col[col_ptr[i]];
    }
    
    for( i = 0; i < local_m; i++ ){
        Node_rc_t *x = prob->rX[i];
        while( x->c_index != -1 ){
            int idx = x->c_index; // without - 1
            X_col[col_ptr[idx]].r_index = i; // starts from 0, (1)
            X_col[col_ptr[idx]].c_index = idx;
            X_col[col_ptr[idx]].value = x->value;
            col_ptr[idx]++;
            x++;
        }
    }

    for( i = 0; i < local_n; i++ ){
        X_col[col_ptr[i]].r_index = -1;
        X_col[col_ptr[i]].c_index = -1;
    }

    free(col_ptr);

}

void read_scaling_problem( struct Problem *prob ){

	int 		max_index;
	char 		*label, *idx, *val;
	size_t 		entries;
	Node_rc_t 	*X_tmp;
	Node_rc_t 	**tmp_rX;

	FILE *fp = fopen( param.prob_name, "r" );

    if( fp == NULL ){
        fprintf( stdout, "Can't open input file %s\n", param.prob_name );
        exit(1);
    }

	max_line_size = 1024;
	problem.num_samples = 0;
	entries = 0;

	line = Malloc( char, max_line_size );

	max_index = -1;

	while( get_line(fp) != NULL ){
	
		label = strtok( line, " \t\n");
        /* features */
        while(1){

            idx = strtok( NULL, ":" );
            val = strtok( NULL, " \t" );

            if( val == NULL || idx == NULL )
                break;

            if( max_index < atoi(idx) - 1  ){
            	max_index = atoi(idx) - 1;
            }

            entries++;      
        }
        problem.num_samples++; 		
	}
	rewind( fp );

	problem.num_features = max_index + 1;

	problem.y = Malloc( double, problem.num_samples );

	X_tmp = Malloc( Node_rc_t, entries + problem.num_samples );

	tmp_rX = Malloc( Node_rc_t *, problem.num_samples );

	fprintf(stdout, "Dataset Description:\n" );
	fprintf(stdout, "Number of samples : %d\n", problem.num_samples );
	fprintf(stdout, "Number of features: %d\n", problem.num_features );

	mean_column = Malloc( double, problem.num_features );
	min_column = Malloc( double, problem.num_features );
	max_column = Malloc( double, problem.num_features );

	int max_entries = problem.num_features * problem.num_samples;

	if( param.scaling == 1 ){
		entries = max_entries;
	}else if( param.scaling == 2 ){
		entries = 100;
	}

	fprintf( stdout, "Collecting statistic...\n" );

	int j = 0;
	int i;

	for ( i = 0; i < problem.num_samples; i++ ){
	
		get_line(fp);

		tmp_rX[i] = &X_tmp[j];
		/* strtok: The C library function char *strtok(char *str, const char *delim) breaks 
		  string str into a series of tokens using the delimiter delim. */
		label = strtok( line, " \t\n");
		
		if( label == NULL )
			continue;

		problem.y[i] = atof(label);

		while(1){
			/* The first call to strtok() returns a pointer to the first token in the string 
			 pointed to by s1. Subsequent calls to strtok() must pass a NULL pointer as the 
			 first argument, in order to get the next token in the string. */
			idx = strtok( NULL, ":" );
			val = strtok( NULL, " \t" );

			if( val == NULL || idx == NULL )
				break;

			int index = atoi(idx) - 1;
			double value = atof(val);

			if( min_column[index] > value ){
				min_column[index] = value;
			}	

			if( max_column[index] < value ){
				max_column[index] = value;
			}

			mean_column[index] += value;

			X_tmp[j].c_index = index; /* c_index is index of feature */
			X_tmp[j].r_index = i; /* r_index is index of sample */
			X_tmp[j].value = value;

			j++;
	
		}	

		X_tmp[j].c_index = -1;
		X_tmp[j++].r_index = -1;

	}

	rewind( fp );

	for( i = 0; i < problem.num_features; i++ ){
		mean_column[i] /= (double) problem.num_samples;
	}

	if( param.scaling == 1 ){
		X = Malloc( Node_rc_t, max_entries + problem.num_samples );
		j = 0;
		for( i = 0; i < problem.num_samples; i++ ){
			int fea_idx = 0;
			problem.rX[i] = &X[j];
			for( fea_idx = 0; fea_idx < problem.num_features; fea_idx++ ){
				X[j+fea_idx].value = ( X[j+fea_idx].value - mean_column[fea_idx] );
			}

			while( tmp_rX[i]->c_index != -1 ){
				tmp_rX[i]++;
			}
			
		}
	}

	/* select sparse or dense data format for scaling type 2 */
	if( param.scaling == 2 ){
		int count2 = 0;
		for( i = 0; i < problem.num_samples; i++ ){
			while( tmp_rX[i]->c_index != -1 ){
				if( tmp_rX[i]->value == min_column[i] ){
					count2++;
				}
				tmp_rX[i]++;
			}
		}
		fprintf(stdout, "Scaling 2 needs to be done\n" );
		fprintf(stdout, "count2: %d\n", count2);
	}

	fclose( fp );
}


/* Output function based on print level*/
void output(){

	FILE *fo;

	/* if print level is zero, print nothing */
	if( param.print_level == 0 ){
	}
	/* if print level is one, print basic result */
	else if( param.print_level == 1 ){
		if( output_farsa.term == 1 ){
			fprintf( stdout, "Optimal solution has been found.\n" );
			fprintf( stdout, "Objective function value: %f\n", F );
			fprintf( stdout, "Iteration: %d\n", iter);
			fprintf( stdout, "Error: %f\n", MAX( norm_beta, norm_phi ) );
			fprintf( stdout, "Target tolerance: %f\n", ttol );
		}else if( output_farsa.term == 0 ){
			fprintf( stdout, "Maximum iteration has been reached.\n" );
		}
		fprintf( stdout, "Run time: %f\n", output_farsa.run_time );
	}
	/* if print level is one, print detailed result on final iteration */
	else if( param.print_level == 2 ){

		/* display summary result in terminal */
		if( output_farsa.term == 1 ){
			fprintf( stdout, "Optimal solution has been found.\n" );
			fprintf( stdout, "Objective function value: %f\n", F );
			fprintf( stdout, "Iteration: %d\n", iter);
			fprintf( stdout, "Error: %2.2e\n", MAX( norm_beta, norm_phi ) );
			fprintf( stdout, "Target tolerance: %2.2e\n", ttol );
		}else if( output_farsa.term == 0 ){
			fprintf( stdout, "Maximum iteration has been reached.\n" );
		}
		fprintf( stdout, "Run time: %f seconds\n", output_farsa.run_time );
		
		/* output results into output file */
		/* if output_file is not set, set it as a temporary file */
		if( param.output_file == NULL ){
			strncpy( param.output_file, "temp_output.txt", sizeof(param.output_file)-1);
		}else{
			int len;
			len = (int) strlen( param.output_file );	
			if( len == 0 ){
				strncpy( param.output_file, "temp_output.txt", sizeof(param.output_file)-1);
			}	
		}
		fo = fopen( param.output_file, "w" );
		if( fo == NULL ){
	        fprintf( stdout, "Can't open output file %s\n", param.output_file );
	        exit(1);		
		}		
		fprintf( fo, "FaRSA output\n\n\n" );

		fprintf( fo, "objective function type: %d\n", param.objective_function_type );
		if( param.objective_function_type != 0 ){
			fprintf( fo, "data file:               %s\n", param.prob_name );
			fprintf( fo, "data format:             %s\n", param.data_format );
			fprintf( fo, "number of features:      %d\n", problem.num_features );
			fprintf( fo, "number of samples :      %d\n", problem.num_samples );		
		}else{
			fprintf( fo, "n:                       %d\n", n );
		}	

		fprintf( fo, "lambda:                  %e\n", lambda );
		fprintf( fo, "target tolerance:        %e\n\n", output_farsa.ttol );

		fprintf( fo, "%s\n", double_hline ); 

		fprintf( fo, "Optimization solver performance information.\n" );

		fprintf( fo, "%s\n", hline ); 

		if( output_farsa.term == 1 ){
			fprintf( fo, "termination status:      Optimal solution has been found.\n" );
		}else if( output_farsa.term == 0 ){
			fprintf( fo, "termination status:      Maximum iteration has been found.\n" );
		}


		fprintf( fo, "initial error:           %e\n", MAX( output_farsa.norm_beta_initial, output_farsa.norm_phi_initial ) );
		fprintf( fo, "final error:             %e\n", output_farsa.error );
		fprintf( fo, "final l2 norm of beta:   %e\n", output_farsa.norm_beta_final );
		fprintf( fo, "final l2 norm of phi:    %e\n", output_farsa.norm_phi_final );
		fprintf( fo, "initial func value:      %f\n", output_farsa.F_initial );
		fprintf( fo, "final func value:        %f\n", output_farsa.F_final );
		fprintf( fo, "function eval:           %d\n", output_farsa.f_eval );
		fprintf( fo, "gradient eval:           %d\n", output_farsa.grad_eval );
		fprintf( fo, "hessVec  eval:           %d\n", output_farsa.hessVec_eval );
		fprintf( fo, "iteration:               %d\n", output_farsa.iter );
		fprintf( fo, "beta iteration:          %d\n", output_farsa.beta_iter );
		fprintf( fo, "phi iteration:           %d\n", output_farsa.phi_iter );
		fprintf( fo, "zero perc in solution:   %f\n", output_farsa.zero_perc );
		fprintf( fo, "\n%s\n", hline );

		fprintf( fo, "run time:                %f seconds\n", output_farsa.run_time );
		fprintf( fo, "%s\n", hline );

		fclose( fo );

	}
	/* if print level is 3, print result in column for run all*/
	else if( param.print_level == 3 ){
		/* display summary result in terminal */	

		fprintf( stdout, "problem, function_eval, gradient_eval, hessVec_eval, zero_perc, " );
		fprintf( stdout, "iter, beta_iter, phi_iter, run_time\n" );
		fprintf( stdout, "%s, %d, %d, %d, %f, ", param.name, output_farsa.f_eval, output_farsa.grad_eval, output_farsa.hessVec_eval, output_farsa.zero_perc );	
		fprintf( stdout, "%d, %d, %d, %f\n", output_farsa.iter, output_farsa.beta_iter, output_farsa.phi_iter, output_farsa.run_time );

	}
	/* if print level is 4, print result of each iteration, since each iteration result has been outputted */
	else if( param.print_level == 4 ){

		/* start final results output */
		fo = fopen( param.output_file, "ab" );
		
		/* Print solver result */
		fprintf( fo, "\n" );
		fprintf( fo, "%s\n", hline );
		fprintf( fo, "Final  result\n" );
		fprintf( fo, "=============\n" );
		if ( output_farsa.term == 0 ) {
			fprintf( fo, "  EXIT: Iteration limit reached\n" );
		} else if ( output_farsa.term == 1 ) {
			fprintf( fo, "  EXIT: Optimal solution found\n" );
		}

		fprintf( fo, "\n" );

		/* Print iterate quantities */
		fprintf( fo, "Final values\n" );
		fprintf( fo, "============\n" );
		fprintf( fo, "  Objective function......................... : %+e\n", output_farsa.f_final );
		fprintf( fo, "  One norm of the solution................... : %+e\n", output_farsa.x_l1_norm_final ) ;
		fprintf( fo, "  Objective function plus one norm........... : %+e\n", output_farsa.F_final );
		fprintf( fo, "  Target tolerance........................... : %+e\n", output_farsa.ttol );
		fprintf( fo, "  Optimality error........................... : %+e\n", MAX( output_farsa.norm_beta_final, output_farsa.norm_phi_final ) );
		fprintf( fo, "  Zero percentage............................ : %+e\n", output_farsa.zero_perc );
		fprintf( fo, "  Status..................................... : %d \n", output_farsa.status );
		fprintf( fo, "\n" );

		/* Print counters */
		fprintf( fo, "Final counters\n" );
		fprintf( fo, "==============\n" );
		fprintf( fo, "  Iterations................................. : %8d  \n", output_farsa.iter );
		fprintf( fo, "  Beta iterations............................ : %8d  \n", output_farsa.beta_iter );
		fprintf( fo, "  Phi iterations............................. : %8d  \n", output_farsa.phi_iter );
		fprintf( fo, "  Function evaluations....................... : %8d  \n", output_farsa.f_eval );
		fprintf( fo, "  Gradient evaluations....................... : %8d  \n", output_farsa.grad_eval );
		fprintf( fo, "  Hessian vector products.................... : %8d  \n", output_farsa.hessVec_eval );		
   		fprintf( fo, "  CPU seconds................................ : %8.4e\n", output_farsa.run_time );


	}
	/* if print level is 5, print result of each iteration, since each iteration result has been outputted */
	else if( param.print_level == 5 ){

		/* display summary result in terminal */
		if( output_farsa.term == 1 ){
			fprintf( stdout, "Optimal solution has been found.\n" );
			fprintf( stdout, "Objective function value: %f\n", F );
			fprintf( stdout, "Iteration: %d\n", iter);
			fprintf( stdout, "Error: %2.2e\n", MAX( norm_beta, norm_phi ) );
			fprintf( stdout, "Target tolerance: %2.2e\n", ttol );
		}else if( output_farsa.term == 0 ){
			fprintf( stdout, "Maximum iteration has been reached.\n" );
		}
		fprintf( stdout, "Run time: %f seconds\n", output_farsa.run_time );

		/* start final results output */
		fo = fopen( param.output_file, "ab" );
		
		/* Print solver result */
		fprintf( fo, "\n" );
		fprintf( fo, "%s\n", hline );
		fprintf( fo, "Final  result\n" );
		fprintf( fo, "=============\n" );
		if ( output_farsa.term == 0 ) {
			fprintf( fo, "  EXIT: Iteration limit reached\n" );
		} else if ( output_farsa.term == 1 ) {
			fprintf( fo, "  EXIT: Optimal solution found\n" );
		}

		fprintf( fo, "\n" );

		/* Print iterate quantities */
		fprintf( fo, "Final values\n" );
		fprintf( fo, "============\n" );
		fprintf( fo, "  Objective function......................... : %+e\n", output_farsa.f_final );
		fprintf( fo, "  One norm of the solution................... : %+e\n", output_farsa.x_l1_norm_final ) ;
		fprintf( fo, "  Objective function plus one norm........... : %+e\n", output_farsa.F_final );
		fprintf( fo, "  Target tolerance........................... : %+e\n", output_farsa.ttol );
		fprintf( fo, "  Optimality error........................... : %+e\n", MAX( output_farsa.norm_beta_final, output_farsa.norm_phi_final ) );
		fprintf( fo, "  Zero percentage............................ : %+e\n", output_farsa.zero_perc );
		fprintf( fo, "  Status..................................... : %d \n", output_farsa.status );
		fprintf( fo, "\n" );

		/* Print counters */
		fprintf( fo, "Final counters\n" );
		fprintf( fo, "==============\n" );
		fprintf( fo, "  Iterations................................. : %8d  \n", output_farsa.iter );
		fprintf( fo, "  Beta iterations............................ : %8d  \n", output_farsa.beta_iter );
		fprintf( fo, "  Phi iterations............................. : %8d  \n", output_farsa.phi_iter );
		fprintf( fo, "  Function evaluations....................... : %8d  \n", output_farsa.f_eval );
		fprintf( fo, "  Gradient evaluations....................... : %8d  \n", output_farsa.grad_eval );
		fprintf( fo, "  Hessian vector products.................... : %8d  \n", output_farsa.hessVec_eval );		
   		fprintf( fo, "  CPU seconds................................ : %8.4e\n", output_farsa.run_time );


	}

}


/* Output of traveling */
void output_problem_data() {

	FILE *fo;
	// fprintf( stdout, "output_problem_data\n" );
	fo = fopen( param.output_file, "ab" );

	fprintf( fo, "Problem Attributes\n" );
	fprintf( fo, "==================\n" );
	fprintf( fo, "  Name....................................... : %s  \n", param.prob_name );
	fprintf( fo, "  Number of variables........................ : %10d  \n", problem.num_features );
	fprintf( fo, "  Number of samples.......................... : %10d  \n", problem.num_samples );
	fprintf( fo, "  Value for lambda........................... : %10.4e\n", lambda );
	fprintf( fo, "\n" );
	fclose(fo);

}

/* Output parameters */
void output_param() {

	FILE *fo;

	fo = fopen( param.output_file, "ab" );

	fprintf( fo, "Control Parameters\n" );
	fprintf( fo, "==================\n" );
	fprintf( fo, " max_iterates............................... : %10d  \n", param.max_iter );
	fprintf( fo, " max_time................................... : %10f  \n", param.max_time);
	fprintf( fo, " print_level................................ : %10d  \n", param.print_level );
	fprintf( fo, " print_file................................. : %s    \n", param.output_file );
	fprintf( fo, " print_every................................ : %10d  \n", param.print_every );
	fprintf( fo, " Gamma...................................... : %10.4e\n", param.Gamma );
	fprintf( fo, " eta_r...................................... : %10.4e\n", param.eta_r );
	fprintf( fo, " eta........................................ : %10.4e\n", param.eta );
	fprintf( fo, " xi......................................... : %10.4e\n", param.xi );
	fprintf( fo, " phiFrac.................................... : %10.4e\n", param.phiFrac );
	fprintf( fo, " betaFrac................................... : %10.4e\n", param.betaFrac );
	fprintf( fo, " sub_solver................................. : %10d  \n", param.sub_solver );
	fprintf( fo, " maxCG_iterates............................. : %10d  \n", param.max_CG_iter );
	fprintf( fo, " maxCDqp_iterates........................... : %10d  \n", param.max_CD_iter );
	fprintf( fo, " CrossOverTol............................... : %10.4e\n", param.crossOverTol );
	fprintf( fo, " maxback.................................... : %10d  \n", param.maxback );
	fprintf( fo, " fracViol................................... : %10.4e\n", param.fracViol );
	fprintf( fo, " nVmaxAllow................................. : %10d  \n", param.nVmaxAllow );
	fprintf( fo, " term_type.................................. : %10d  \n", param.termType );
	fprintf( fo, " absolute tolerance......................... : %10.4e\n", param.tol_absolute );



	fclose( fo );

}

void output_header() {

	FILE *fo;
	
	fo = fopen( param.output_file, "ab" );

	fprintf( fo, "%s\n", hline );
	fprintf( fo, "%s\n", column_titles );
	fprintf( fo, "%s\n", hline );
	
	fclose( fo );
}

void free_all(){

	if( I_phi != 0 )
		free(I_phi);
	if( I_z != 0 )
		free(I_z);
	if( I_n != 0 )
		free(I_n);
	if( I_p != 0 )
		free(I_p);
	if( beta != 0 )
		free(beta);	
	if( phi != 0 )
		free(phi);
	if( grad_f != 0 )
		free(grad_f);
	if( grad_F != 0 )
		free(grad_F);
	if( grad_f != 0 )
		free(grad_f);
	if( x != 0 )
		free(x);
	if( y != 0 )
		free(y);
	if( d != 0 )
		free(d);
	if( x_linesearch != 0 )
		free(x_linesearch);
	if( step != 0 )
		free(step);
	if( X != 0 )
		free(X);
	if( X_col != 0 )
		free(X_col);
	if( X_row_S != 0 )
		free(X_row_S);
	if( col_ptr != 0 )
		free(col_ptr);
	if( col_ptr_head != 0 )
		free(col_ptr_head);	
	if( mean_column != 0 )
		free(mean_column);	
	if( min_column != 0 )
		free(min_column);
	if( max_column != 0 )
		free(max_column);		

}

void swap( int *a, int *b ){
	
	int tmp;

	tmp = *a;
	*a 	= *b;
	*b 	= tmp;

}

void copy( double *vector1, double *vector2, int n ) {

	int i;

	for ( i = 0; i < n; i++ ) {
		vector1[i] = vector2[i];
	}
}
