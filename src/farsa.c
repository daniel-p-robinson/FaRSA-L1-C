/**** 
	FaRSA

	version: 2.0

	Authors:

	Tianyi Chen, Frank E. Curtis, Daniel P. Robinson
	
	April, 2018
*/
#include "farsa.h"

Output_FaRSA_t farsa( int argc, char **argv, Input_FaRSA_t *input_farsa )
{

	/* cpu time */
	clock_t 	begin_cpu;
	clock_t 	end_cpu;
	double 		cpu_time;

	/* get arguments from parser */
	parse_command_line( argc, argv );

	/* load setting from profile */
	load_setting_profile();

	/* initialize rest parameters */
	/* load dataset is involved here. */
	initialize_rest_param( *input_farsa );

	/* if loss function has optimized version, then directly call optimized routine
	  if not, run the generic routine. */	
	begin_cpu 	= clock();

	/* training */
	switch ( param.objective_function_type ) 
	{
		case LOGISTIC_LOSS:
			if( param.print_level >= 1 && param.print_level != 3 && param.print_level != 4 ){
				fprintf( stdout, "Logistic loss plus l1 regularizer...\n" );
			}
			logistic_loss( &problem );
			break;
		case ELASTIC_NET:
			if( param.print_level >= 1 && param.print_level != 3 && param.print_level != 4 ){
				fprintf( stdout, "Least square loss plus l1 regularizer...\n" );
			}
			elastic_net_loss( &problem );
			break;			
		default:

			if( param.print_level >= 1 ){
				fprintf( stdout, "objective function type is invalid \n" );
			}
			exit(0);		
	}

	end_cpu 				= clock();
	cpu_time 				= (double)(end_cpu - begin_cpu) / CLOCKS_PER_SEC;

	output_farsa.run_time 	= cpu_time;
	/* output based on print level */
	output();

	/* run test routine */
	if ( param.do_test == TRUE ) 
	{

		Test_output_t 	test_result;
		Problem_t 		problem_test;

		/* read test data file */
		read_sparse_problem( &problem_test, FALSE );

		/* do testing */
		test_result 	= test( &problem_test );

		/* output_test */	
		output_test( test_result );		

	}	

	if ( param.print_final_x == TRUE ) {
		output_final_x();
	}

	return output_farsa;

}


/************** Read parses and settings ************/
void parse_command_line( int argc, char **argv )
{

	int i;
	/* set default values */
	param.max_iter 		= 1000;
	param.print_level 	= 2;
	param.Gamma 		= 1.0;
	param.eta_r 		= 1E-1;
	param.eta 			= 1E-2;
	param.xi 			= 0.5;
	param.tol_absolute 	= 1E-6;
	param.tol_relative 	= 1E-6;
	param.betaFrac 		= 0.8;
	param.phiFrac 		= 1.0;
	param.tryCG 		= 0;
	param.tryCD 		= 1;
	param.tryCrossOver 	= 0;
	param.max_CG_iter 	= MAX_VALUE;
	param.max_CD_iter 	= 1000;
	param.crossOverTol 	= 1E-1;
	param.TRradiusPhi 	= 1E2;
	param.TRradiusBeta 	= 1E-2;
	param.maxback 		= 100;
	param.fracViol 		= 0.25;
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
	param.do_test 		= FALSE;
	param.print_final_x = FALSE;


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
void load_setting_profile()
{

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
			strncpy( param.train_file, (char *)value, sizeof(param.train_file)-1 );
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
		} else if ( !strcmp( attribute, "do_test" ) ) {
			param.do_test = atoi( value );
			continue;
		} else if ( !strcmp( attribute, "test_data_file" ) ) {
			int len;
			len = (int) strlen( (char *)value );
			if( len > 1024 ){
				fprintf(stdout, "test data file path is too long\n" );
				exit(-1);
			}
			strncpy( param.test_file, (char *)value, sizeof(param.test_file)-1 );
			continue;
		} else if ( !strcmp( attribute, "print_final_x" ) ) {
			param.print_final_x = atoi( value );
			continue;
		}
	}

	fclose(fp);
}

/* get line function */
char* get_line( FILE *input )
{

	int len;

	if( fgets( line, max_line_size, input ) == NULL )
	{
		return NULL;
	}

	/* if current max_line_size doesn't find the \n, then enlarge 
	 the max_line_size */
	while( strrchr( line, '\n' ) == NULL )
	{
		max_line_size *= 2;
		line = (char *) realloc( line, max_line_size );
		len = (int) strlen(line);
		if( fgets( line+len, max_line_size-len, input ) == NULL )
		{
			break;
		}
	}
	return line;

}

/* initialize parameters after loading profile */
void initialize_rest_param( struct Input_FaRSA input_farsa  )
{

	/* initialize parameters */
	/* objective function is not personalized, then read problem */
	if( param.objective_function_type != 0 )
	{
		/* read libsvm format dataset */
		read_sparse_problem( &problem, TRUE );

		m 			= problem.num_samples;
		n 			= problem.num_features;
		y 			= Malloc( double, m );
		m_d 		= (double) m; /* double type of m */
		n_d 		= (double) n;
		lambda 		= param.lambda < 0 ? param.lambdaCoeff / (double) problem.num_samples : param.lambda;
		
	}
	else if( param.objective_function_type == 0 )
	{

		if( input_farsa.n == 0 )
		{
			fprintf(stdout, "Current n is 0. Please set n. n is the number of variables.\n" );
			exit(1);
		}
		lambda 		= param.lambda;
		n 			= input_farsa.n;
		m 			= 1;
		m_d 		= (double) m;

	}
	else
	{
		
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
	y_B 			= Malloc( double, n );
	dirDer 			= 0.0;
	F_old 			= 0.0;  /* objective function value for last step */
	F 				= 0.0; /* objective function value on current step */
	step 			= Malloc( double, n ); /* difference between two steps for trust region update */
	norm_step 		= 0.0;
	maxback 		= param.maxback;
	xi 				= param.xi;
	sameOrthant 	= TRUE;

	row_ptr_sp 		= Malloc( int, m + 1 );
	help_hv_logis 	= Malloc( double, m );
}


/* read libsvm format sparse dataset, following by a similar function in liblinear */
void read_sparse_problem( Problem_t *prob, int is_train )
{

	int i;
	int j;
	int max_index;
	char *label, *idx, *val;
	size_t entries;
	FILE *fp;

	if ( is_train == TRUE ) 
	{
		fp = fopen( param.train_file, "r" );
	} 
	else 
	{
		fp = fopen( param.test_file, "r" );
	}

	
    if( fp == NULL )
    {
        fprintf( stdout, "Can't open input file %s\n", param.test_file );
        // exit(0) behave like return 0 in main() function, exit(1) behave like return 1 
        exit(1);
    }

	max_line_size = 1024;
	prob->num_samples = 0;
	entries = 0;

	
	line = Malloc( char, max_line_size );

	while( get_line(fp) != NULL )
	{

		label = strtok( line, " \t\n");

        // features
        while(TRUE)
        {

            idx = strtok( NULL, ":" );
            val = strtok( NULL, " \t" );

            if( val == NULL || idx == NULL )
            {
                break;
            }

            entries++;      
        }
        prob->num_samples++; 		
	}
	rewind( fp );

	prob->y 		= Malloc( double, prob->num_samples );
	prob->rX 		= Malloc( Node_t *, prob->num_samples ); 
	prob->nnz_rows 	= Malloc( int, prob->num_samples );
	prob->nnz 		= entries;

	/* sort by samples */
	X_row = Malloc( Node_t, entries + prob->num_samples );

	max_index = -1;

	j = 0;

	for ( i = 0; i < prob->num_samples; i++ )
	{
	
		get_line(fp);

		prob->rX[i] 		= &X_row[j];
		prob->nnz_rows[i] 	= 0;

		/* strtok: The C library function char *strtok(char *str, const char *delim) breaks 
		 string str into a series of tokens using the delimiter delim. */
		label = strtok( line, " \t\n");
		
		if( label == NULL )
		{
			continue;
		}

		prob->y[i] = atof(label);

		while(TRUE)
		{
			/* The first call to strtok() returns a pointer to the first token in the string 
			 pointed to by s1. Subsequent calls to strtok() must pass a NULL pointer as the 
			 first argument, in order to get the next token in the string. */
			idx = strtok( NULL, ":" );
			val = strtok( NULL, " \t" );

			if( val == NULL || idx == NULL )
			{
				break;
			}

			X_row[j].index 		= atoi(idx) - 1; /* index for node in X_row is index of feature */
			X_row[j].value 		= atof(val);
            
            /* if the current index of feature is larger than maximum one, update the maximum index of feature */
            if( max_index < X_row[j].index )
            {
                max_index = X_row[j].index;
            }
			j++;

			prob->nnz_rows[i]++;
		}	
		X_row[j].index 		= -1;
		j++;
	}
	prob->num_features 	= max_index + 1;

	if ( is_train == TRUE ) { 
		transpose( prob );
	} 

	if( param.print_level >= 1 && param.print_level != 3 && param.print_level != 4 ){
		fprintf(stdout, "Dataset Description:\n" );
		fprintf(stdout, "Number of samples : %d\n", prob->num_samples );
		fprintf(stdout, "Number of features: %d\n", prob->num_features );		
	}

	fclose(fp);
}

/* transpose row-wise of X, rX into column-wise of X, cX */
void transpose( struct Problem *prob )
{

    int i;
    
    int local_m;
    int local_n;

    local_m = prob->num_samples;
    local_n = prob->num_features;
    
    int nnz;
    nnz = 0;

    col_ptr = Malloc( int, local_n + 1 );


	if ( col_ptr == NULL) 
	{
		fprintf( stdout, "Memory issue, please wait for a few seconds and run again! \n" );
        exit(-1);
    }

    
    col_ptr_head 	= Malloc( int, local_n + 1 );
    prob->cX 		= Malloc( Node_t *, prob->num_features );
    prob->nnz_cols 	= Malloc( int, local_n );

    for( i = 0; i < local_n + 1; i++ )
    {
        col_ptr[i] 			= 0;
        col_ptr_head[i] 	= 0;
        prob->nnz_cols[i] 	= 0;
    }

    for( i = 0; i < local_m; i++ ){
        Node_t *x = prob->rX[i];
        while( x->index != -1 ){
            nnz++;
            col_ptr[x->index+1]++;
            col_ptr_head[x->index+1]++;
            prob->nnz_cols[x->index]++;
            x++;
        }
    }

    
    for( i = 1; i < local_n + 1; i++ ){
        col_ptr[i] += col_ptr[i-1] + 1;
        col_ptr_head[i] += col_ptr_head[i-1] + 1;
    }

    X_col = Malloc( Node_t, nnz + local_n );

    for( i = 0; i < local_n; i++ ){
        prob->cX[i] = &X_col[col_ptr[i]];
    }
    
    for( i = 0; i < local_m; i++ ){
        Node_t *x = prob->rX[i];
        while( x->index != -1 ){
            int idx = x->index; // the feature index
            X_col[col_ptr[idx]].index = i; /* index of sample */
            X_col[col_ptr[idx]].value = x->value;
            col_ptr[idx]++;
            x++;
        }
    }

    for( i = 0; i < local_n; i++ ){
    	X_col[col_ptr[i]].index 	= -1;
    }

    free(col_ptr);

    /* for v2.0, free X_row */
    free(X_row);
    free(prob->rX);
    // free(prob->nnz_rows);
}


/* set beta and phi */
void setBetaPhi() 
{	

	int 			i;
	int 			*p_phi;
	int 			*p_beta;

	p_phi 			= I_phi;
	p_beta 			= I_beta;

	norm_beta 		= 0.0;
	norm_phi 		= 0.0;

	nnz 			= 0;

	for( i = 0; i < n; i++ )
	{
		beta[i] = 0.0;
		phi[i] = 0.0;
	}


	for ( i = 0; i < n; i++ )
	{

		if( x[i] == 0 )
		{
			I_z[i] = 1;
			I_p[i] = 0;
			I_n[i] = 0;
			if ( grad_f[i] < -lambda )
			{
				beta[i] = grad_f[i] + lambda;			
			}
			else if( grad_f[i] > lambda )
			{
				beta[i] = grad_f[i] - lambda;
			}
			/* use l2 norm at first */
			norm_beta += beta[i] * beta[i];
			phi[i] = 0.0;

			if ( beta[i] != 0.0 ) 
			{
				nI_beta++;
				*p_beta = i;
				p_beta++;
			}

		}
		else if( x[i] > 0 )
		{
			phi[i] = grad_f[i] + lambda;
			if( phi[i] > 0 )
			{
				phi[i] = MIN( phi[i], MAX( x[i], grad_f[i] - lambda ) );
			}
			norm_phi += phi[i] * phi[i];
			I_p[i] = 1;
			I_z[i] = 0;
			I_n[i] = 0;
			grad_F[i] = grad_f[i] + lambda; 
			if( phi[i] != 0.0 )
			{
				nI_phi++;
				*p_phi = i;
				p_phi++;
			}
			
			nnz++;

		}
		else
		{
			phi[i] = grad_f[i] - lambda;
			if( phi[i] < 0 )
			{
				phi[i] = MAX( phi[i], MIN( x[i], grad_f[i] + lambda ) );
			}
			norm_phi += phi[i] * phi[i];
			I_n[i] = 1;
			I_z[i] = 0;
			I_p[i] = 0;
			grad_F[i] = grad_f[i] - lambda;
			if( phi[i] != 0.0 )
			{
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
}


/* select beta */
void select_beta() 
{
	if ( param.betaFrac >= 1.0 ) 
	{
		Beta 		= beta;
		norm_Beta 	= norm_beta;
		return ;
	} 
	else 
	{	

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
		for ( i = 0; i < n; i++ ) 
		{
			Beta[i] 	= 0.0;
		}

		norm_Beta 	= 0.0;
		for ( i = 0; i < nI_beta_selected; i++ ) 
		{
			Beta[I_beta_selected[i]]	= beta[I_beta_selected[i]];
			norm_Beta += Beta[I_beta_selected[i]] * Beta[I_beta_selected[i]];
		}

		norm_Beta 	= sqrt( norm_Beta );

	}
}

/* select phi */
void select_phi() 
{

	if ( param.phiFrac >= 1.0 ) 
	{
		I_phi_selected 		= I_phi;
		nI_phi_selected 	= nI_phi;
		return ;

	} 
	else 
	{

		int i;
		phi_I 			= Malloc( double, nI_phi );
		I_phi_selected 	= Malloc( int, nI_phi );

		for ( i = 0; i < nI_phi; i++ ) 
		{
			phi_I[i] 			= phi[I_phi[i]];
			I_phi_selected[i]	= I_phi[i];
		}
		nI_phi_selected = MIN( nnz, MAX( ceil( (double) nnz * param.phiFrac ), 5000 ) );
		quick_select( phi_I, nI_phi_selected - 1, 0, nI_phi - 1, I_phi_selected );
		quick_sort( I_phi_selected, 0, nI_phi_selected - 1 );
	
	}

}

/* quick select largest K from vector */
double quick_select( double *vector, int K, int lo, int hi, int *indexes ) 
{

	int i;
	while ( lo <= hi ) 
	{
		i = partition_descend( vector, lo, hi, indexes );
		if ( i < K ) 
		{
			lo = i + 1;
		} 
		else if ( i > K ) 
		{
			hi = i - 1;
		} 
		else 
		{
			return vector[K];
		}
	}
	return vector[K];

}

int partition_descend( double *vector, int lo, int hi, int *indexes ) 
{
	int 	i, j;
	double 	pivot;
	double 	temp;
	int  	temp_idx;

	pivot 	= fabs(vector[lo]);
	i 		= lo;
	j 		= hi;

	while( TRUE ) 
	{

		while ( fabs(vector[i]) >= pivot && i <= hi ) 
		{
			i++;
		}

		while ( fabs(vector[j]) < pivot ) 
		{
			j--;
		}

		if ( i >= j ) 
		{
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
	indexes[lo] = temp_idx;
	
	return j;

}

/* quick sort with duplicate keys by three-way partition */
void quick_sort( int* array, int lo, int hi ) 
{
	
	int median_idx;
	int pivot_left_idx;
	int pivot_right_idx;
	int pivot;
	int i;
	if ( hi <= lo ) 
	{	
		return;
	}
	median_idx 		= lo + ( hi - lo ) / 2;
	pivot_left_idx 	= lo;
	pivot_right_idx = hi;
	swap( &array[lo], &array[median_idx] );

	pivot 			= array[lo];
	i 				= lo;

	while ( i <= pivot_right_idx ) 
	{
		if ( array[i] < pivot ) 
		{
			swap( &array[pivot_left_idx], &array[i] );
			pivot_left_idx++;
			i++;
		} 
		else if ( array[i] > pivot ) 
		{
			swap( &array[i], &array[pivot_right_idx] );
			pivot_right_idx--;
		} 
		else 
		{
			i++;
		}		
	}

	quick_sort( array, lo, pivot_left_idx - 1 );
	quick_sort( array, pivot_right_idx + 1, hi );

}

/***********  Logistic Regression Loss Section Start *********/
/**
 *	The routine contains 
 */
void logistic_loss( Problem_t *prob )
{
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
	Input_CD_t 		input_cd;
	Output_sp_t 	output_cp;

	/* help_vector for logistic loss is sigmoid vector */
	diagonal 		= Malloc( double, m );
	sigmoid 		= Malloc( double, m );
	ysigmoid 		= Malloc( double, m );
	expterm 		= Malloc( double, m );
	wTx 			= Malloc( double, m );


	/* Set labels and Initialize some expterm ysigmoid */
	for ( i = 0; i < m; i++ ) 
	{
		if ( problem.y[i] > 0 ) 
		{
			y[i] 			= 1.0;
		} 
		else 
		{
			y[i] 			= -1.0;
		}
		expterm[i] 			= 1.0;
		sigmoid[i] 			= 0.5;
		ysigmoid[i] 		= y[i] * sigmoid[i];
		diagonal[i] 		= 0.25;
		wTx[i] 				= 0.0;
		help_hv_logis[i] 	= 0.0;
	}

	/* Continue initialize gradient */
	for( i = 0 ; i < n; i++ )
	{
		x[i] 		= 0.0;
		grad_f[i] 	= -1.0 / m_d * dot_ds( ysigmoid, prob->cX[i] );
	}	

	grad_eval++;

	func_output 	= logistic_func( x, expterm );
	F 				= func_output.F_value;
	f_eval++;

	/* Initially set beta phi */
	setBetaPhi();

	norm_phi0 		= norm_phi;
	norm_beta0 		= norm_beta;

	/* set target termination tolerance */
	if( param.termType == 2 )
	{
		ttol = MAX( tol_absolute, tol_relative * MAX( norm_beta0, norm_phi0 ) );
	}
	else if( param.termType == 1 )
	{
		ttol = tol_absolute;
	}

	/* set a few attributes in output of farsa */
	output_farsa.x_initial 			= x;
	output_farsa.norm_beta_initial 	= norm_beta0;
	output_farsa.norm_phi_initial 	= norm_phi0;
	output_farsa.ttol 				= ttol;		

	/* print to output file if print level >= 4 */
	if ( param.print_level >= 4 ) 
	{

		fo = fopen( param.output_file, "w" );
		if ( fo == NULL ) 
		{
			fprintf( stdout, "Can not open %s\n", param.output_file );
			exit(0);
		}

		fprintf( fo, "%s\n", hline );
		fprintf( fo, "|                  A Fast Reduced Space Algorithm (FaRSA) for Sparse Logistic Regression v. 2.0                     |\n" );
		fprintf( fo, "%s\n", double_hline );

		fclose( fo );

		output_problem_data();
		
		output_param();
	}	

	/* Main loop for FaRSA for logistic loss */
	while( TRUE ) 
	{	
		/* print based on print level */
		if ( param.print_level >= 4 ) 
		{

			fo = fopen( param.output_file, "ab" );
			
			if ( iter % param.print_every == 0 ) 
			{

				fprintf( fo, "%s\n", hline_1 );
				fprintf( fo, "%s\n", column_titles );
				fprintf( fo, "%s\n", hline_1 );			
			}
			
			fprintf( fo, "%5d %9.2e %8.2e %9.5e %8.2e %8.2e %5d |", iter, func_output.f_value, func_output.x_l1_norm, func_output.F_value, norm_phi, norm_beta, nnz );
			
			fclose( fo );

		}

		/* termination condition */
		/* terminate if error is smaller than the tolerance */
		if( MAX( norm_beta, norm_phi ) <= ttol )
		{

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
			output_farsa.n 				 = n;
			break;
		}

		/* terminate if current iteration exceed maximum allowed iteration */
		if( iter >= param.max_iter )
		{
			
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
			output_farsa.n 				 = n;
			break;
		
		}

		/* if norm of beta is less than the norm of phi, then do phi step */
		if( norm_beta <= Gamma * norm_phi )
		{

			sprintf( type_iteration, "%s", "phi" );
			
			phi_iter++;
			
			iter_type = 1;

			select_phi();

			/* solve sub-problem based on corresponding type */
			switch ( param.sub_solver ) 
			{

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
					input_cg.d_full 	 = d;
					input_cg.hessVecProd = &logistic_hessVecProd;

					/* set up data for phi-step */
					logistic_setXF();

					output_cp 			 = CGsolver( input_cg );

					free(problem.cXS);

					break;

				/* Compute the search direction using CD on the quadratic model. */
				case USE_CD_QP:

					break;
				default:

					break;	

			} /* subproblem solver ends */

			/* Perform a line search */
			nV 					= output_cp.nV;
			j 					= 0;
			alpha 				= 1.0;
			nDifforthant 		= 0;
			dirDer 				= output_cp.dir_der;
			sp_residual 		= output_cp.res;
			sp_target_residual 	= output_cp.res_target;
			norm_d 				= output_cp.norm_d;
			nVmax 				= output_cp.nVmax;
			suf_descent_1 		= param.eta * dirDer; 
			F_old 				= F;
			proj_output 		= project( x, d, alpha, x_linesearch );
			nDifforthant 		= proj_output.nV;
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

			/* check whether current iterate close enough to x^**/
			if( nV == 0 && norm_beta < tol_absolute )
			{
				rsType = 2;
			}
			else
			{
				rsType = 1;
			}

			/* while loop for line search */
			while( TRUE ) 
			{

				norm_step 	= 0.0;

				if( !sameOrthant && F < F_old )
				{
					/* set x, and grad_f */
					for( i = 0; i < n; i++ )
					{
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						x[i] = x_linesearch[i];
						grad_f[i] = -1.0 / m_d * dot_ds( ysigmoid, problem.cX[i] );
					}
					grad_eval++;
					norm_step = sqrt( norm_step );
					sprintf( linesearch_flag, "%s", "orthant" );
					break;
				}


				if( sameOrthant && F - F_old <=  suf_descent_1 * alpha )
				{
					
					if ( sameOrthant_prev || alpha > MIN( 1E-4, ttol ) ) 
					{
						
						for( i = 0; i < n; i++ )
						{
							step[i] = x_linesearch[i] - x[i];
							norm_step += step[i] * step[i];
							x[i] = x_linesearch[i];
							grad_f[i] = -1.0 / m_d * dot_ds( ysigmoid, problem.cX[i] );
						}
					
						grad_eval++;
						norm_step = sqrt( norm_step );
						sprintf( linesearch_flag, "%s", "descent" );
						break;
					
					} 
					else 
					{

						alpha_B 		= MAX_VALUE;

						for ( i = 0; i < n; i++ ) 
						{

							/* x[i] / d[i] < 0.0 */
							if ( d[i] < 0.0 && x[i] > 0.0 ) 
							{
								alpha_B = MIN( alpha_B, -x[i] / d[i] );
							}

							if ( d[i] > 0.0 && x[i] < 0.0 ) 
							{
								alpha_B = MIN( alpha_B, -x[i] / d[i] );
							}

						}


						for ( i = 0; i < n; i++ ) 
						{							
							y_B[i] 		= x[i] + alpha_B * d[i];
						}

						logistic_setExpTerm( y_B, wTx, sigmoid, ysigmoid, expterm );
						func_output 		= logistic_func( y_B, expterm );
						F_B 				= func_output.F_value;
						f_eval++;

						if ( F_B - F_old <= alpha_B * suf_descent_1 ) 
						{
						
							x_linesearch = y_B;
							F = F_B;
							/* set x, and grad_f */
							for( i = 0; i < n; i++ )
							{
								step[i] = x_linesearch[i] - x[i];
								norm_step += step[i] * step[i];
								x[i] = x_linesearch[i];
								grad_f[i] = -1.0 / m_d * dot_ds( ysigmoid, problem.cX[i] );
							}
							grad_eval++;
							norm_step = sqrt( norm_step );
							break;

						} 
						else 
						{
							
							logistic_setExpTerm( x_linesearch, wTx, sigmoid, ysigmoid, expterm );
							/* set x, and grad_f */
							for( i = 0; i < n; i++ )
							{
								step[i] = x_linesearch[i] - x[i];
								norm_step += step[i] * step[i];
								x[i] = x_linesearch[i];
								grad_f[i] = -1.0 / m_d * dot_ds( ysigmoid, problem.cX[i] );
							}
							grad_eval++;
							norm_step = sqrt( norm_step );
							break;
						
						}

					}
				
				}

				if ( j > maxback )
				{
					/* set x, and grad_f */
					for( i = 0; i < n; i++ )
					{
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						
						x[i] = x_linesearch[i];
						grad_f[i] = -1.0 / m_d * dot_ds( ysigmoid, problem.cX[i] );
					
					}

					grad_eval++;
					norm_step = sqrt( norm_step );
					break;
				}		

				alpha 				*= xi;
				proj_output 		= project( x, d, alpha, x_linesearch );
				nDifforthant 		= proj_output.nV;
				//x_linesearch 		= proj_output.project_vector;
				sameOrthant_prev 	= sameOrthant;
				sameOrthant 		= proj_output.same_sign;
				logistic_setExpTerm( x_linesearch, wTx, sigmoid, ysigmoid, expterm );
				func_output 		= logistic_func( x_linesearch, expterm );
				F 					= func_output.F_value;
				f 					= func_output.f_value;
				x_l1_norm 			= func_output.x_l1_norm;
				f_eval++;
				j++;

			} /* while loop for line search ends */

		}
		else
		{
			beta_iter++;
			select_beta();

			iter_type 			= 0;
			TRbeta_norm_beta 	= -( TRradiusBeta / norm_Beta );

			sprintf( type_iteration, "%s", "beta" );
			sub_iters 			= 1;

			for ( i = 0; i < n; i++ )
			{
				d[i] = TRbeta_norm_beta * Beta[i];
				x_linesearch[i] = x[i] + d[i];
				if( d[i] > 0 )
				{
					grad_F[i] = grad_f[i] + lambda;
				}
				else if( d[i] < 0 )
				{
					grad_F[i] = grad_f[i] - lambda;
				}
				else
				{
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

			/* perform line search */
			while( TRUE )
			{

				norm_step 		= 0.0;
				
				if( F - F_old <= alpha * suf_descent_1 )
				{
					/* set x, and grad_f */
					for( i = 0; i < n; i++ )
					{
						step[i] 	= x_linesearch[i] - x[i];
						norm_step 	+= step[i] * step[i];
						x[i] 		= x_linesearch[i];
						grad_f[i] 	= -1.0 / m_d * dot_ds( ysigmoid, problem.cX[i] );
					}
					grad_eval++;
					norm_step 		= sqrt( norm_step );
					sprintf( linesearch_flag, "%s", "descent" );
					break;
				}	

				if ( j > maxback )
				{
					/* set x, and grad_f */
					for( i = 0; i < n; i++ )
					{
						step[i] 	= x_linesearch[i] - x[i];
						norm_step 	+= step[i] * step[i];
						x[i] 		= x_linesearch[i];
						grad_f[i] 	= -1.0 / m_d * dot_ds( ysigmoid, problem.cX[i] );
					}
					grad_eval++;
					norm_step 		= sqrt( norm_step );
					break;
				}
				j++;
				alpha 		*= xi;	

				/* set new trial step */
				for( i = 0 ; i < n; i++ )
				{
					x_linesearch[i] = x[i] + alpha * d[i];
				}

				logistic_setExpTerm( x_linesearch, wTx, sigmoid, ysigmoid, expterm );
				func_output = logistic_func( x_linesearch, expterm );
				F 			= func_output.F_value;
				f 			= func_output.f_value;
				x_l1_norm 	= func_output.x_l1_norm;
				f_eval++;

			}		


		} /* end beta step */

		/* output based on different iteration type */
		if ( param.print_level >= 4 ) 
		{

			fo = fopen( param.output_file, "ab" );

			if ( !strcmp( type_iteration, "phi" ) ) 
			{
				// fprintf(stdout, "outputting \n" );
				fprintf( fo, " %4s %5d %7.1e %7s %6d %5d %8.2e %8.2e %5d %5d %8.2e |", type_iteration, nI_phi, TRradiusPhi, output_cp.sub_prob_flag, rsType, sub_iters, sp_residual, sp_target_residual, output_cp.nV, nVmax, norm_d );
				fprintf( fo, " %3d %8.2e %7s %5d |\n", j, alpha, linesearch_flag, nDifforthant );
				// %4d %4d %1.3f %1.3f | %1.3f %1.3f %1.3f\n",  j, alpha,  );
			} 
			else if ( !strcmp( type_iteration, "beta" ) ) 
			{
				fprintf( fo, " %4s %5d %7.1e ------- %6d %5d -------- -------- ----- ----- %8.2e |", type_iteration, nI_phi, TRradiusBeta, rsType, sub_iters, norm_d );
				fprintf( fo, " %3d %8.2e %7s %5d |\n", j, alpha, linesearch_flag, 0 );
			}

			fclose( fo );
		}

		for ( i = 0; i < n; i++ ) 
		{
			d[i] = 0.0;
		}

		nI_phi 	= 0;
		nI_beta	= 0;
		
		/* Set trust-region radius for next iteration. */
		if( iter_type == 1 )
		{
			TRradiusPhi 	= MAX( 1E-1, MIN( 1E3, 10 * norm_step ) );
		}
		else
		{
			TRradiusBeta 	= MAX( 1E-5, MIN( 1.0, norm_step ) );
		}
		
		norm_step = 0.0;

		setBetaPhi();

		iter++;

	} /* end of main loop of farsa logistic cost function */


}

/* set terms for logistic cost */
void logistic_setExpTerm( double *x, double *wTx, double *sigmoid, double *ysigmoid, double *expterm )
{
	
	int i;
   	Node_t *tmp_x_col;

   	for ( i = 0; i < n; i++ )
   	{
   		tmp_x_col 	= problem.cX[i];
   		/* the index of tmp_x_col is the sample index */
   		while( tmp_x_col->index != -1 )
   		{
   			wTx[tmp_x_col->index] 	+= x[i] * tmp_x_col->value;
   			tmp_x_col++;
   		}
   	}

   	for( i = 0; i < m; i++ ) 
   	{
    	expterm[i] 	= exp( -1.0 * y[i] * wTx[i] );
    	sigmoid[i] 	 = expterm[i] / ( 1.0 + expterm[i] );
    	ysigmoid[i] = y[i] * sigmoid[i]; 
    	diagonal[i] = sigmoid[i] * ( 1.0 - sigmoid[i] );
		/* remember to set wTx as 0 here, so next time, no need to call another for loop */		
    	wTx[i] 		= 0.0; 
    }	   
}

/* compute logisitic loss function given x and expterm */
Func_output_t logistic_func( double *x, double *expterm )
{
	
	Func_output_t 	func_output;
 	int 			i;
   
	func_output.f_value = 0.0;
	for ( i = 0; i < m; i++ ) 
	{
		func_output.f_value += log( 1.0 + expterm[i] );
	}

	func_output.f_value 	/= (double) m;
	func_output.x_l1_norm 	= l1_n5( x, n );
	func_output.F_value 	= func_output.f_value + lambda * func_output.x_l1_norm;

	return func_output;

}


/* set data for phi-step subproblem */
void logistic_setXF()
{
	int 	i;

	/* allocate column format data for reduced space problem */
	problem.cXS = Malloc( Node_t *, nI_phi_selected );

	for( i = 0; i < nI_phi_selected; i++ )
	{
		problem.cXS[i] = &X_col[col_ptr_head[I_phi_selected[i]]];
	}	
}

/* logistic loss hess vector product */
void logistic_hessVecProd( double *v, double *hv ){

	int 	i;
	Node_t  *tmp_x_col;


	for ( i = 0; i < m; i++ )
	{
		help_hv_logis[i] 	= 0.0;
	}

	for ( i = 0; i < nI_phi_selected; i++ )
	{
		tmp_x_col 	= problem.cXS[i];
		/* the index of tmp_x_col is the index of sample */
		while ( tmp_x_col->index != -1 )
		{	
			help_hv_logis[tmp_x_col->index] += diagonal[tmp_x_col->index] * tmp_x_col->value * v[i];	
			tmp_x_col++;
		}
	}

	for( i = 0; i < nI_phi_selected; i++ ) {
		hv[i] = 1.0 / m_d * dot_ds( help_hv_logis, problem.cXS[i] ) + 1E-8 * v[i];				
	}

}

/* optimized routine for least-square loss */
void elastic_net_loss( struct Problem *prob )
{
	int 			i;
	int 			j;
	int 			nDifforthant;
	int 			nV;
	int 			nVmax;
	int 			n_threshold;

	double 			suf_descent_1;
	double 			initial_F;
	double 			norm_d;
	double 			TRbeta_norm_beta;
	double 			tmp;
	double 			sp_residual;
	double 			sp_target_residual;
	double 			T;
	double 			cost_save;
	double 			cost_not_save;
	FILE 			*fo;
	char 			linesearch_flag[40];

	lambda_l2 		= 0.0;

	Node_t 			*curr_node_i;
	Node_t 			*curr_node_j;
	Func_output_t	func_output;
	Proj_output_t 	proj_output;
	Input_CG_t 		input_cg;
	Input_CD_t 		input_cd;
	Output_sp_t 	output_cp;

	ATy 			= Malloc( double, n ); 
	Ax 				= Malloc( double, m );
	AI_x 			= Malloc( double, m );
	y 				= prob->y;	
	norm_y_sq 		= 0.0;
	n_threshold 	= 10000;
	T 				= 10.0;

	/* determine whether precomputing A^TA or not */
	if ( n >= n_threshold )
	{
		saving_flag = NOT_SAVE;
	} 
	else
	{
		saving_flag = MAY_SAVE;

		// double nnz_d;
		// nnz_d 		= (double) prob->nnz;

		// cost_save 	= T * n_d * n_d - 1.0 / 3.0 * T * n_d + ( m_d + 1.0 ) / 3.0 * nnz_d;

		// cost_not_save 	= 4.0 / 3.0 * T * nnz_d;

		// for ( i = 0; i < m; i++ )
		// {
		// 	cost_save 	+= ( m_d + 1.0 ) / 6.0 * (double)( MAX( 0, prob->nnz_cols[i] - 1 ) );
		// }
		// fprintf(stderr, "%f\n", cost_save);
		if ( 2 * prob->nnz > ( n * n ) )
		{
			saving_flag = MUST_SAVE;
		}
	}

	/* for test */
	// saving_flag 	= NOT_SAVE;
	saving_flag 	= MUST_SAVE;

	/* if saving flag is must_save, then pre-compute A^TA */
	if ( saving_flag == MUST_SAVE )
	{
		ATA 		= Malloc( double *, n );
		for ( i = 0; i < n; i++ )
		{	
			ATA[i] 	= Malloc( double, n );
		}
		/* for each row */
		for ( i = 0; i < n; i++ )
		{
			curr_node_i 		= prob->cX[i];
			/* for each column */
			for ( j = 0; j < n; j++ )
			{	
				curr_node_j 	= prob->cX[j];
				ATA[i][j] 		= dot_ss( curr_node_i, curr_node_j );		
			}
		}		
	}

	/* Continue initialize gradient */
	for( i = 0 ; i < n; i++ ){
		x[i] 		= 0.0;
		ATy[i] 		= dot_ds( y, prob->cX[i] );
		grad_f[i] 	= -ATy[i] + 2.0 * lambda_l2 * x[i];
	}	

	grad_eval++;

	if ( saving_flag == MUST_SAVE )
	{
		for ( i = 0; i < m; i++ )
		{
			Ax[i] 		= 0.0;
			norm_y_sq 	+= y[i] * y[i];	
		}		
	}
	else
	{
		for ( i = 0; i < m; i++ )
		{
			Ax[i] 		= 0.0;
			norm_y_sq 	+= y[i] * y[i];	
			// AI_x[i] 	= 0.0;
		}		
	}


	func_output 	= least_square_func( x );
	F 				= func_output.F_value;
	f_eval++;

	/* Initially set beta phi */
	setBetaPhi();

	norm_phi0 		= norm_phi;
	norm_beta0 		= norm_beta;

	/* set target termination tolerance */
	if ( param.termType == 2 )
	{
		ttol = MAX( tol_absolute, tol_relative * MAX( norm_beta0, norm_phi0 ) );
	}
	else if ( param.termType == 1 )
	{
		ttol = tol_absolute;
	}

	/* set a few attributes in output of farsa */
	output_farsa.x_initial 			= x;
	output_farsa.norm_beta_initial 	= norm_beta0;
	output_farsa.norm_phi_initial 	= norm_phi0;
	output_farsa.ttol 				= ttol;		

	/* print to output file if print level >= 4 */
	if ( param.print_level >= 4 ) 
	{

		fo = fopen( param.output_file, "w" );
		if ( fo == NULL ) 
		{
			fprintf( stdout, "Can not open %s\n", param.output_file );
			exit(0);
		}

		fprintf( fo, "%s\n", hline );
		fprintf( fo, "|                  A Fast Reduced Space Algorithm (FaRSA) for Least Square Loss v. 2.0                     |\n" );
		fprintf( fo, "%s\n", double_hline );

		fclose( fo );

		output_problem_data();
		
		output_param();
	}	

	/* Main loop for FaRSA for logistic loss */
	while( TRUE ) 
	{	
		/* print based on print level */
		if ( param.print_level >= 4 ) {

			fo = fopen( param.output_file, "ab" );
			
			if ( iter % param.print_every == 0 ) 
			{

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
			output_farsa.n 				 = n;
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
			output_farsa.n 				 = n;
			break;
		
		}

		/* if norm of beta is less than the norm of phi, then do phi step */
		if( norm_beta <= Gamma * norm_phi )
		{

			sprintf( type_iteration, "%s", "phi" );
			
			phi_iter++;
			
			iter_type = 1;

			select_phi();

			/* solve sub-problem based on corresponding type */
			switch ( param.sub_solver ) 
			{

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
					input_cg.d_full 	 = d;
					input_cg.hessVecProd = &least_square_hessVecProd;

					/* set up data for phi-step */
					if ( saving_flag == NOT_SAVE )
					{
						least_square_setXF();				
					}

					output_cp 			 = CGsolver( input_cg );

					// free(problem.cXS);
					break;

				/* Compute the search direction using CD on the quadratic model. */
				case USE_CD_QP:

					break;
				default:

					break;	

			} /* subproblem solver ends */

			/* Perform a line search */
			nV 					= output_cp.nV;
			j 					= 0;
			alpha 				= 1.0;
			nDifforthant 		= 0;
			dirDer 				= output_cp.dir_der;
			sp_residual 		= output_cp.res;
			sp_target_residual 	= output_cp.res_target;
			norm_d 				= output_cp.norm_d;
			nVmax 				= output_cp.nVmax;
			suf_descent_1 		= param.eta * dirDer; 
			F_old 				= F;

			proj_output 		= project( x, d, alpha, x_linesearch );
			nDifforthant 		= proj_output.nV;
			sameOrthant 		= proj_output.same_sign;
			sameOrthant_prev	= TRUE;

			least_sqaure_setAx( x_linesearch );
			func_output 		= least_square_func( x_linesearch );
			F 					= func_output.F_value;
			f 					= func_output.f_value;
			x_l1_norm 			= func_output.x_l1_norm;
			sub_iters 			= output_cp.iter;
			sprintf( linesearch_flag, "%s", output_cp.sub_prob_flag );
			f_eval++;

			/* check whether current iterate close enough to x^**/
			if( nV == 0 && norm_beta < tol_absolute )
			{
				rsType = 2;
			}
			else
			{
				rsType = 1;
			}

			/* while loop for line search */
			while( TRUE ) 
			{

				norm_step 	= 0.0;

				if( !sameOrthant && F < F_old )
				{
					/* set x, and grad_f */
					for( i = 0; i < n; i++ )
					{
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						x[i] = x_linesearch[i];
						if ( saving_flag == MUST_SAVE )
						{
							grad_f[i] 	= dot_n5( ATA[i], x, n ) - ATy[i] + 2.0 * lambda_l2 * x[i];
						} 
						else
						{
							grad_f[i] 	= dot_ds( Ax, prob->cX[i] ) - ATy[i] + 2.0 * lambda_l2 * x[i];
						}
					}
					grad_eval++;
					norm_step = sqrt( norm_step );
					sprintf( linesearch_flag, "%s", "orthant" );
					break;
				}


				if( sameOrthant && F - F_old <=  suf_descent_1 * alpha )
				{
					
					if ( sameOrthant_prev || alpha > MIN( 1E-4, ttol ) ) 
					{
						
						for( i = 0; i < n; i++ )
						{
							step[i] = x_linesearch[i] - x[i];
							norm_step += step[i] * step[i];
							x[i] = x_linesearch[i];
							if ( saving_flag == MUST_SAVE )
							{
								grad_f[i] 	= dot_n5( ATA[i], x, n ) - ATy[i] + 2.0 * lambda_l2 * x[i];
							} 
							else
							{
								grad_f[i] 	= dot_ds( Ax, prob->cX[i] ) - ATy[i] + 2.0 * lambda_l2 * x[i];
							}
						}
					
						grad_eval++;
						norm_step = sqrt( norm_step );
						sprintf( linesearch_flag, "%s", "descent" );
						break;
					
					} 
					else 
					{

						alpha_B 		= MAX_VALUE;

						for ( i = 0; i < n; i++ ) 
						{

							/* x[i] / d[i] < 0.0 */
							if ( d[i] < 0.0 && x[i] > 0.0 ) 
							{
								alpha_B = MIN( alpha_B, -x[i] / d[i] );
							}

							if ( d[i] > 0.0 && x[i] < 0.0 ) 
							{
								alpha_B = MIN( alpha_B, -x[i] / d[i] );
							}

						}


						for ( i = 0; i < n; i++ ) 
						{							
							y_B[i] 		= x[i] + alpha_B * d[i];
						}

						least_sqaure_setAx( y_B );
						func_output 		= least_square_func( y_B );
						F_B 				= func_output.F_value;
						f_eval++;

						if ( F_B - F_old <= alpha_B * suf_descent_1 ) 
						{
						
							x_linesearch = y_B;
							F = F_B;
							/* set x, and grad_f */
							for( i = 0; i < n; i++ ){
								step[i] = x_linesearch[i] - x[i];
								norm_step += step[i] * step[i];
								x[i] = x_linesearch[i];
								if ( saving_flag == MUST_SAVE )
								{
									grad_f[i] 	= dot_n5( ATA[i], x, n ) - ATy[i] + 2.0 * lambda_l2 * x[i];
								} 
								else
								{
									grad_f[i] 	= dot_ds( Ax, prob->cX[i] ) - ATy[i] + 2.0 * lambda_l2 * x[i];
								}
							}
							grad_eval++;
							norm_step = sqrt( norm_step );
							break;

						} 
						else 
						{
							least_sqaure_setAx( x_linesearch );
							/* set x, and grad_f */
							for( i = 0; i < n; i++ )
							{
								step[i] = x_linesearch[i] - x[i];
								norm_step += step[i] * step[i];
								x[i] = x_linesearch[i];
								if ( saving_flag == MUST_SAVE )
								{
									grad_f[i] 	= dot_n5( ATA[i], x, n ) - ATy[i] + 2.0 * lambda_l2 * x[i];
								} 
								else
								{
									grad_f[i] 	= dot_ds( Ax, prob->cX[i] ) - ATy[i] + 2.0 * lambda_l2 * x[i];
								}
							}
							grad_eval++;
							norm_step = sqrt( norm_step );
							break;
						
						}

					}
				
				}

				if ( j > maxback )
				{
					/* set x, and grad_f */
					for( i = 0; i < n; i++ )
					{
						step[i] = x_linesearch[i] - x[i];
						norm_step += step[i] * step[i];
						
						x[i] = x_linesearch[i];

						if ( saving_flag == MUST_SAVE )
						{
							grad_f[i] 	= dot_n5( ATA[i], x, n ) - ATy[i] + 2.0 * lambda_l2 * x[i];
						} 
						else
						{
							grad_f[i] 	= dot_ds( Ax, prob->cX[i] ) - ATy[i] + 2.0 * lambda_l2 * x[i];
						}
					
					}

					grad_eval++;
					norm_step = sqrt( norm_step );
					break;
				}		

				alpha 				*= xi;
				proj_output 		= project( x, d, alpha, x_linesearch );
				nDifforthant 		= proj_output.nV;
				sameOrthant_prev 	= sameOrthant;
				sameOrthant 		= proj_output.same_sign;
				least_sqaure_setAx( x_linesearch );
				func_output 		= least_square_func( x_linesearch );
				F 					= func_output.F_value;
				f 					= func_output.f_value;
				x_l1_norm 			= func_output.x_l1_norm;
				f_eval++;
				j++;

			} /* while loop for line search ends */

		}
		else
		{
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

			least_sqaure_setAx( x_linesearch );
			func_output 		= least_square_func( x_linesearch );
			F 					= func_output.F_value;
			f 					= func_output.f_value;
			x_l1_norm 			= func_output.x_l1_norm;	
			f_eval++;

			suf_descent_1 		= param.eta * dirDer; 
			tmp 				= 0.0; /* used for update grad_f */

			/* perform line search */
			while( TRUE ){

				norm_step 		= 0.0;
				
				if( F - F_old <= alpha * suf_descent_1 ){
					/* set x, and grad_f */
					for( i = 0; i < n; i++ ){
						step[i] 	= x_linesearch[i] - x[i];
						norm_step 	+= step[i] * step[i];
						x[i] 		= x_linesearch[i];

						if ( saving_flag == MUST_SAVE )
						{
							grad_f[i] 	= dot_n5( ATA[i], x, n ) - ATy[i] + 2.0 * lambda_l2 * x[i];
						} 
						else
						{
							grad_f[i] 	= dot_ds( Ax, prob->cX[i] ) - ATy[i] + 2.0 * lambda_l2 * x[i];
						}

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

						if ( saving_flag == MUST_SAVE )
						{
							grad_f[i] 	= dot_n5( ATA[i], x, n ) - ATy[i] + 2.0 * lambda_l2 * x[i];
						} 
						else
						{							
							grad_f[i] 	= dot_ds( Ax, prob->cX[i] ) - ATy[i] + 2.0 * lambda_l2 * x[i];
						}

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

				least_sqaure_setAx( x_linesearch );
				func_output = least_square_func( x_linesearch );
				F 			= func_output.F_value;
				f 			= func_output.f_value;
				x_l1_norm 	= func_output.x_l1_norm;
				f_eval++;

			}		


		} /* end beta step */

		/* output based on different iteration type */
		if ( param.print_level >= 4 ) 
		{

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

		for ( i = 0; i < n; i++ ) 
		{
			d[i] = 0.0;
		}

		nI_phi 	= 0;
		nI_beta	= 0;
		
		/* Set trust-region radius for next iteration. */
		if( iter_type == 1 )
		{
			TRradiusPhi 	= MAX( 1E-1, MIN( 1E3, 10 * norm_step ) );
		}
		else
		{
			TRradiusBeta 	= MAX( 1E-5, MIN( 1.0, norm_step ) );
		}
		
		norm_step = 0.0;

		setBetaPhi();

		iter++;

	} /* end of main loop for least-square loss + l1 */

}

/* calculate least-square loss function value */
Func_output_t least_square_func( double *x )
{
	Func_output_t 	func_output;
   	
   	int 	i;

	func_output.f_value = 0.0;

	for ( i = 0; i < m; i++ )
	{
		func_output.f_value += ( Ax[i] - y[i] ) *  ( Ax[i] - y[i] );
	}
	func_output.f_value 	*= 0.5;
	func_output.x_l1_norm 	= l1_n5( x, n );
	func_output.F_value 	= func_output.f_value + lambda_l2 * dot_n5( x, x, n ) + lambda * func_output.x_l1_norm;

	return func_output;
}

/* update Ax */
void least_sqaure_setAx( double *x )
{
	int i;
	Node_t *curr_node;

	for ( i = 0; i < m; i++ )
	{
		Ax[i] 	= 0.0;
	}

	for ( i = 0; i < n; i++ )
	{
		curr_node 	= problem.cX[i];

		while ( curr_node->index != -1 )
		{
			Ax[curr_node->index] += x[i] * curr_node->value;
			curr_node++;
		}
	}
}

/* least square loss set data for phi-step subproblem */
void least_square_setXF()
{
	int 	i;

	/* allocate column format data for reduced space problem */
	problem.cXS = Malloc( Node_t *, nI_phi_selected );

	for( i = 0; i < nI_phi_selected; i++ )
	{
		problem.cXS[i] = &X_col[col_ptr_head[I_phi_selected[i]]];
	}	
}

/* least square loss hess vector product */
void least_square_hessVecProd( double *v, double *hv )
{
	int 	i;
	int 	j;
	int 	idx_1;
	int 	idx_2;


	if ( saving_flag == MUST_SAVE )
	{
		for ( i = 0; i < nI_phi_selected; i++ )
		{
			idx_1 	= I_phi_selected[i];
			hv[i] 	= 2.0 * lambda_l2 * v[i];

			for ( j = 0; j < nI_phi_selected; j++ )
			{
				idx_2 	= I_phi_selected[j];
				hv[i] 	+= ATA[idx_1][idx_2] * v[j];  
			} 
		}
	} 
	else 
	{
		Node_t 	*curr_node;		
		/* compute AI_X at first */
		for ( i = 0; i < m; i++ )
		{
			AI_x[i] 	= 0.0;
		}

		for ( i = 0; i < nI_phi_selected; i++ )
		{
			curr_node 	=  problem.cXS[i];
			while( curr_node->index != -1 )
			{
				AI_x[curr_node->index] 	+= curr_node->value * v[i];
				curr_node++;
			}
		}

		for ( i = 0; i < nI_phi_selected; i++ )
		{
			hv[i] 	= dot_ds( AI_x, problem.cXS[i] ) + 2.0 * lambda_l2 * v[i];
		}

	}

}


/* calculate sparsity in x */
double calculate_sparsity( double *x, int n )
{

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
/* Operator 1: the inner product of one vector and a sparse vector with the same dimension */
double dot_ds( const double *v, const Node_t *x )
{
   
   double 	vTx;
   
   vTx = 0.0;
   
   while( x->index != -1 ){
      vTx += v[x->index] * x->value;
      x++;
   }
   
   return vTx;
}

/* Operator 2: the inner product of two sparse vectors */
double dot_ss( const Node_t *v_1, const Node_t *v_2 )
{
	double v_1Tv_2;

	v_1Tv_2 	= 0.0;

	while( v_1->index != -1 && v_2->index != -1 )
	{
		if ( v_1->index == v_2->index )
		{
			v_1Tv_2 += v_1->value * v_2->value;
			v_1++;
			v_2++;
		}
		else if ( v_1->index > v_2->index )
		{
			v_2++;
		}
		else
		{
			v_1++;
		}
	}

	return v_1Tv_2;
}


/* n5 operators start */
double dot_n5( double *v1, double *v2, int n )
{
   
   int 		i, n5;
   double 	result;
   
   result 	= 0.0;

   if( n <= 0 ) 
   {
   		return result;
   }

   n5 		= n % 5;

   for( i = 0; i < n5; i++ )
   {
      result += v1[i] * v2[i];
   }

   for( ; i < n; i += 5 )
   {
      result += v1[i]*v2[i] + v1[i+1]*v2[i+1] + v1[i+2]*v2[i+2] + v1[i+3]*v2[i+3] + v1[i+4]*v2[i+4];
   }

   return result;
}

double l1_n5( double *v, int n )
{
   int 		i, n5;
   double 	result;

   result 	= 0.0;

   if( n <= 0 ) 
   {
   		return result;
   } 

   n5 		= n % 5;

   for( i = 0; i < n5; i++ )
   {
      result += fabs(v[i]);
   }

   for( ; i < n; i += 5 )
   {
      result += fabs(v[i]) + fabs(v[i+1]) + fabs(v[i+2]) + fabs(v[i+3]) + fabs(v[i+4]);
   }
   return result;
}


/* swap two integers */
void swap( int *a, int *b )
{
	
	int tmp;

	tmp = *a;
	*a 	= *b;
	*b 	= tmp;

}

/* CG reduced-space solver */
Output_sp_t CGsolver( Input_CG_t input_cg ){

	Output_sp_t output_cg;
	
	int 		i;
	int 		nV;
	int 		iter_cg;
	int 		*I;
	int 		nI;
	int 		local_n;
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

	nV 			= 0;
	iter_cg 	= 0;
	I 			= input_cg.I;
	nI 			= input_cg.nI;
	local_n 	= input_cg.n;
	nVmaxAllow 	= input_cg.nVmaxAllow;
	max_CG_iter = input_cg.maxCG_iter;
	max_loop 	= MIN( max_CG_iter,nI );
	x 			= input_cg.x;
	grad_F 		= input_cg.grad_F;
	TRradius 	= input_cg.TRradius;
	d_full 		= input_cg.d_full;
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
	

	for( i = 0; i < local_n; i++ ){
		d_full[i] = 0.0;
	}


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
		input_cg.hessVecProd( p, Hp );
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
		for( i = 0; i < local_n; i++ ){
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
	output_cg.nV 			= nV;
	output_cg.iter 			= iter_cg;
	output_cg.res 			= norm_r;
	output_cg.norm_d 		= norm_d;
	output_cg.res_target 	= res_target;
	output_cg.nVmax 		= nVmaxAllow;
	sprintf( output_cg.sub_prob_flag, "%s", sub_prob_flag );

	free(d_reduced);
	free(r);
	free(p);
	free(Hp);
	return output_cg;	

}



/****************** functions that can be used for any function including special optimized loss function ***/
Proj_output_t project( double *x, double *d, double alpha, double *proj_x )
{

	int 			i;
	int 			same_sign;
	int 			nV;
	Proj_output_t 	proj_output;
	
	same_sign 		= TRUE;
	nV 				= 0;

	for( i = 0; i < n; i++ )
	{
	
		proj_x[i] = x[i] + alpha * d[i];
	
		if( I_p[i] == 1 && proj_x[i] < 0 )
		{
	
			proj_x[i] = 0.0;
			same_sign = FALSE;
			nV++;

		}
		else if( I_n[i] == 1 && proj_x[i] > 0 )
		{
	
			proj_x[i] = 0.0;
			same_sign = FALSE;
			nV++;

		}

	}
	
	proj_output.same_sign 		= same_sign;
	proj_output.nV 				= nV;

	return proj_output;
}

/* Output function based on print level*/
void output(){

	FILE *fo;

	/* if print level is zero, print nothing */
	if( param.print_level == 0 )
	{
		/* pass print nothing */
	}
	/* if print level is 1, print basic result on console */
	else if( param.print_level == 1 )
	{
		if( output_farsa.term == 1 )
		{
			fprintf( stdout, "Optimal solution has been found.\n" );
			fprintf( stdout, "Objective function value: %f\n", F );
			fprintf( stdout, "Iteration: %d\n", iter);
			fprintf( stdout, "Error: %2.2e\n", MAX( norm_beta, norm_phi ) );
			fprintf( stdout, "Target tolerance: %2.2e\n", ttol );
		}
		else if( output_farsa.term == 0 )
		{
			fprintf( stdout, "Maximum iteration has been reached.\n" );
		}
		fprintf( stdout, "Run time: %f seconds\n", output_farsa.run_time );
	}
	/* if print level is 2, print basic result on console and save more detailed information into outputfile */
	else if( param.print_level == 2 )
	{
		/* display summary result in terminal */
		if( output_farsa.term == 1 )
		{
			fprintf( stdout, "Optimal solution has been found.\n" );
			fprintf( stdout, "Objective function value: %f\n", F );
			fprintf( stdout, "Iteration: %d\n", iter);
			fprintf( stdout, "Error: %2.2e\n", MAX( norm_beta, norm_phi ) );
			fprintf( stdout, "Target tolerance: %2.2e\n", ttol );
		}
		else if( output_farsa.term == 0 )
		{
			fprintf( stdout, "Maximum iteration has been reached.\n" );
		}
		fprintf( stdout, "Run time: %f seconds\n", output_farsa.run_time );
		
		/* output results into output file */
		/* if output_file is not set, set it as a temporary file */
		if( param.output_file[0] == '\0' )
		{
			strncpy( param.output_file, "temp_output.txt", sizeof(param.output_file)-1);
		}
		else
		{
			int len;
			len = (int) strlen( param.output_file );	
			if( len == 0 )
			{
				strncpy( param.output_file, "temp_output.txt", sizeof(param.output_file)-1);
			}	
		}
		fo = fopen( param.output_file, "w" );
		if( fo == NULL )
		{
	        fprintf( stdout, "Can't open output file %s\n", param.output_file );
	        exit(1);		
		}		
		fprintf( fo, "FaRSA output\n\n\n" );

		fprintf( fo, "objective function type: %d\n", param.objective_function_type );
		if( param.objective_function_type != 0 )
		{
			fprintf( fo, "data file:               %s\n", param.train_file );
			fprintf( fo, "data format:             %s\n", param.data_format );
			fprintf( fo, "number of features:      %d\n", problem.num_features );
			fprintf( fo, "number of samples :      %d\n", problem.num_samples );		
		}
		else
		{
			fprintf( fo, "n:                       %d\n", n );
		}	

		fprintf( fo, "lambda:                  %e\n", lambda );
		fprintf( fo, "target tolerance:        %e\n\n", output_farsa.ttol );

		fprintf( fo, "%s\n", double_hline ); 

		fprintf( fo, "Optimization solver performance information.\n" );

		fprintf( fo, "%s\n", hline ); 

		if( output_farsa.term == 1 )
		{
			fprintf( fo, "termination status:      Optimal solution has been found.\n" );
		}
		else if( output_farsa.term == 0 )
		{
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
	/* if print level is 3, print final result in column for run all*/
	else if( param.print_level == 3 )
	{
		/* display summary result in console */	

		fprintf( stdout, "problem, function_eval, gradient_eval, hessVec_eval, zero_perc, " );
		fprintf( stdout, "iter, beta_iter, phi_iter, run_time\n" );
		fprintf( stdout, "%s, %d, %d, %d, %f, ", param.name, output_farsa.f_eval, output_farsa.grad_eval, output_farsa.hessVec_eval, output_farsa.zero_perc );	
		fprintf( stdout, "%d, %d, %d, %f\n", output_farsa.iter, output_farsa.beta_iter, output_farsa.phi_iter, output_farsa.run_time );

	}
	/* if print level is 4, print very detailted results into output file, but display nothing on console */
	else if( param.print_level == 4 )
	{

		/* start final results output */
		fo = fopen( param.output_file, "ab" );
		
		/* Print solver result */
		fprintf( fo, "\n" );
		fprintf( fo, "%s\n", hline );
		fprintf( fo, "Final  result\n" );
		fprintf( fo, "=============\n" );
		if ( output_farsa.term == 0 ) 
		{
			fprintf( fo, "  EXIT: Iteration limit reached\n" );
		} 
		else if ( output_farsa.term == 1 ) 
		{
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

   		fclose(fo);
	}
	/* if print level is 5, print very detailted results into output file, and display basic information on console */
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

   		fclose(fo);
	
	}

}

/* Output of traveling */
void output_problem_data() 
{

	FILE *fo;
	// fprintf( stdout, "output_problem_data\n" );
	fo = fopen( param.output_file, "ab" );

	fprintf( fo, "Problem Attributes\n" );
	fprintf( fo, "==================\n" );
	fprintf( fo, "  Name....................................... : %s  \n", param.train_file );
	fprintf( fo, "  Number of variables........................ : %10d  \n", problem.num_features );
	fprintf( fo, "  Number of samples.......................... : %10d  \n", problem.num_samples );
	fprintf( fo, "  Value for lambda........................... : %10.4e\n", lambda );
	fprintf( fo, "\n" );
	fclose(fo);

}

/* Output parameters */
void output_param() 
{

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

void output_header() 
{

	FILE *fo;
	
	fo = fopen( param.output_file, "ab" );

	fprintf( fo, "%s\n", hline );
	fprintf( fo, "%s\n", column_titles );
	fprintf( fo, "%s\n", hline );
	
	fclose( fo );
}


/* test routine */
Test_output_t test( Problem_t *problem_test ) {

	int 			i;
	int 			true_pos_count;
	int 			true_neg_count;
	int 			false_pos_count;
	int 			false_neg_count;
	double 			curr_real_label;
	double 			curr_pred_label;
	double 			accuracy;
	double 			precision;
	double 			recall;
	double 			f1_score;
	double 			*x_final;
	Node_t 			*curr_d;
	Test_output_t 	test_output;

	if ( problem_test->num_features > problem.num_features ) {
		fprintf( stderr, "Warning: Number of features in test file is larger than that of training file, test n: %d, train n: %d\n", problem_test->num_features, problem.num_features );

		x_final 		= Malloc( double, problem_test->num_features );
		copy( x_final, output_farsa.x_final, problem.num_features );		
		/* padding 0 into x_final */
		for ( i = problem.num_features; i < problem_test->num_features; i++ ) {
			x_final[i] 	= 0.0;
		}

		problem_test->num_features 	= problem.num_features;

	} else {
		x_final 		= Malloc( double, problem.num_features );
		copy( x_final, output_farsa.x_final, problem.num_features );
		problem_test->num_features 	= problem.num_features;
	}



	true_pos_count 	= 0;
	true_neg_count 	= 0;
	false_pos_count = 0;
	false_neg_count = 0;

	for ( i = 0; i < problem_test->num_samples; i++ ) {
		
		curr_real_label 	= problem_test->y[i];
		curr_d 				= problem_test->rX[i];

		/* compute p_pos, and p_neg based on different objective function types */
		switch ( param.objective_function_type ) {

			case GENERIC_FUNC:
				break;
			
			case LOGISTIC_LOSS:
				
				curr_pred_label = logistic_predict_one_instance( output_farsa.x_final, curr_d, 0.0 );
				
				if ( curr_real_label > 0 ) {
					curr_real_label = 1;
				} else {
					curr_real_label = -1;
				}
				if ( curr_real_label == 1 && curr_pred_label == 1 ) {
					true_pos_count++;
				} else if ( curr_real_label == -1 && curr_pred_label == -1 ) {
					true_neg_count++;
				} else if ( curr_real_label == 1 && curr_pred_label == -1 ) {
					false_neg_count++;
				} else if ( curr_real_label == -1 && curr_pred_label == 1 ) {
					false_pos_count++;
				}

				break;
			
			case ELASTIC_NET:
			
				break;

			default:

				if( param.print_level >= 1 ){
					fprintf( stdout, "objective function type is invalid \n" );
				}
				exit(0);	
		}

	}

	test_output.accuracy 	= (double)( true_pos_count + true_neg_count ) \
							/ (double)( true_pos_count + true_neg_count + false_pos_count + false_neg_count );
	test_output.precision 	= (double)( true_pos_count ) / (double)( true_pos_count + false_pos_count );
	test_output.recall 		= (double)( true_pos_count ) / (double)( true_pos_count + false_neg_count );
	test_output.f1_score 	= 2.0 * test_output.precision * test_output.recall / ( test_output.precision + test_output.recall );
	test_output.num_samples = problem_test->num_samples;
	test_output.tp 			= true_pos_count;
	test_output.tn  		= true_neg_count;
	test_output.fp 			= false_pos_count;
	test_output.fn 			= false_neg_count;

	return test_output;

}

/* predict label by logisitc regression */
double logistic_predict_one_instance( double *x, const Node_t *d, double bias ) {
	
	double 		p_pos; /* probability to be positive */
	double 		p_neg; /* probability to be negative */
	double 		pred_label;

	p_pos 	= 1.0 / ( 1.0 + exp( -1.0 * ( dot_ds( x, d ) + bias ) ) );
	p_neg 	= 1.0 / ( 1.0 + exp( 1.0 * ( dot_ds( x, d ) + bias ) ) );

	if ( p_pos >= p_neg ) {
		pred_label 	= 1.0;
	} else {
		pred_label 	= -1.0;
	}

	return pred_label;
}

/* output test results */
void output_test( Test_output_t test_result ) {

	FILE *fo;

	/* if print level is zero, print nothing */
	if( param.print_level == 0 ){
	}	
	/* if print level is one, print basic result */
	else if( param.print_level == 1 ){

		fprintf( stdout, "------------------------------\n" );
		fprintf( stdout, "Test Result Summary\n");
		fprintf( stdout, "Number of samples: %d\n", test_result.num_samples );
		fprintf( stdout, "Accuracy:  		 %+e\n", test_result.accuracy );
		fprintf( stdout, "Precision: 		 %+e\n", test_result.precision );
		fprintf( stdout, "Recall: 		     %+e\n", test_result.recall );
		fprintf( stdout, "F1_score:          %+e\n", test_result.f1_score );
	}
	/* if print level is 4, print result of each iteration, since each iteration result has been outputted */
	else if( param.print_level == 4 ){
		/* start final test results output */
		fo = fopen( param.output_file, "ab" );
		
		/* Print solver result */
		fprintf( fo, "\n" );
		fprintf( fo, "%s\n", hline );
		fprintf( fo, "Final test result\n" );
		fprintf( fo, "=============\n" );
		fprintf( fo, "  Number of samples.......................... : %d \n", test_result.num_samples );
		fprintf( fo, "  True positive.............................. : %d \n", test_result.tp );
		fprintf( fo, "  True negative.............................. : %d \n", test_result.tn );
		fprintf( fo, "  False positive............................. : %d \n", test_result.fp );
		fprintf( fo, "  False negative............................. : %d \n", test_result.fn );
		fprintf( fo, "  Accuracy................................... : %+e\n", test_result.accuracy );
		fprintf( fo, "  Precision.................................. : %+e\n", test_result.precision );
		fprintf( fo, "  Recall..................................... : %+e\n", test_result.recall );
		fprintf( fo, "  F1 score................................... : %+e\n", test_result.f1_score );
	
		fclose(fo);
	}
	/* if print level is 5, print result of each iteration, since each iteration result has been outputted */
	else if( param.print_level == 5 ){

		fprintf( stdout, "%s\n", hline );
		fprintf( stdout, "Test Result Summary\n");
		fprintf( stdout, "Number of samples: %d\n", test_result.num_samples );
		fprintf( stdout, "Accuracy:  		 %+e\n", test_result.accuracy );
		fprintf( stdout, "Precision: 		 %+e\n", test_result.precision );
		fprintf( stdout, "Recall: 		     %+e\n", test_result.recall );
		fprintf( stdout, "F1_score:          %+e\n", test_result.f1_score );

		/* start final test results output */
		fo = fopen( param.output_file, "ab" );
		
		/* Print solver result */
		fprintf( fo, "\n" );
		fprintf( fo, "%s\n", hline );
		fprintf( fo, "Final test result\n" );
		fprintf( fo, "=============\n" );
		fprintf( fo, "  Number of samples.......................... : %d \n", test_result.num_samples );
		fprintf( fo, "  True positive.............................. : %d \n", test_result.tp );
		fprintf( fo, "  True negative.............................. : %d \n", test_result.tn );
		fprintf( fo, "  False positive............................. : %d \n", test_result.fp );
		fprintf( fo, "  False negative............................. : %d \n", test_result.fn );
		fprintf( fo, "  Accuracy................................... : %+e\n", test_result.accuracy );
		fprintf( fo, "  Precision.................................. : %+e\n", test_result.precision );
		fprintf( fo, "  Recall..................................... : %+e\n", test_result.recall );
		fprintf( fo, "  F1 score................................... : %+e\n", test_result.f1_score );

		fclose(fo);
	}	
}

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

void print_sparse_vec( const  Node_t *start_node ) {

	while( start_node->index != -1 ){
		fprintf( stdout, "%d:%.2f ", start_node->index, start_node->value );
		start_node++;
	}		
	fprintf(stdout, "\n" );
}

/* output final x */
void output_final_x() {

	FILE *fo;

	fo = fopen( param.output_file, "ab" );

	fprintf( fo, "\n\n" );

	/* Print final x */
	fprintf( fo, "Final x:\n" );
	fprintf( fo, "==============\n" );
	fprintf( fo, "x_final: " );  

	int i;

	for ( i = 0; i < output_farsa.n; i++ ) {
		fprintf( fo, "%e ", output_farsa.x_final[i] );
	}

	fprintf( fo, "\n" );

	fclose(fo);

}

/* copy the elements in vector2 to vector1 */
void copy( double *vector1, double *vector2, int n ) 
{
	int i;
	for ( i = 0; i < n; i++ ) {
		vector1[i] = vector2[i];
	}
}
