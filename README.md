#Fast Reduced Space Algorithm v2.0

FaRSA is an optimized library for solving l1-regularized convex optimization problem:

minimize f(x) + lambda * ||x||_1

## Installation Guide

1. Download source file:

- Clone it from our github repository:

	git clone https://github.com/daniel-p-robinson/FaRSA
	
- Directly download the zip file from our github repository:


2. Install on linux, Mac OSX

Make sure ``gcc" has been installed. Then in FaRSA directory, run the below on terminal:

	make

To test whether installation is successful, run the below command:

	./farsa -t test


## Documentation:

### 1. Quick Start 

In order to successfully call FaRSA to solve aimed problem, users are expected to finish the below steps:

- Configure FaRSA profile

- Prepare readable dataset if necessary

Once the above steps are complete, run the below command in terminal to call FaRSA for solving aimed problem:

	./farsa -p path_to_param_profile

e.g. 
	
	./farsa -p test/FaRSA_t.profile

Results would be displayed in terminal.


### 2. Configure FaRSA Profile

To employ FaRSA on aimed problem correctly, users need to set a FaRSA profile. One demo of profile has been provided in test directory. In profile, 


- objective_function_type: type of f(x), 

	0  generic loss (in test)

	1  logistic loss (well tested)

	2  least square  (in test)

- lambda: value of lambda

	e.g. 0.5 means lambda = 0.5
	
	-1 means lambda is defined by the following lambda_coefficient

- lambda_coefficient: another of define lambda

	lambda = C / num_of_samples

	lambda_coefficient is the C in the above equation

	e.g. 1 means lambda = 1 / num_of_samples

- data_file: the path for data set

- data_format: format of data set, only support libsvm now.

- output_file: results output file path


- print_level: print level

	0  print nothing
	
	1  print basic result on console
	
	2  print basic result on console and save more detailed information into output_file
	
	3  print final result in column for run all
	
	4  print very detailted results into output file, but display nothing on console
	
	5  print very detailted results into output file, and display basic information on console

- tolerance: target tolerance

- term_type: termination type, 

    1  absolute tolerance
    
    2  relative tolerance

- Gamma, eta_r, eta, xi, beta_fraction, phi_fraction, max_CG_iter, max_CD_iter, TRradiusPhi, TRradiusBeta, maxback, frac_viol, nV_max_allow, sub_problem_solver 

	are not recommended to be changed since they are set by parameter tuning. 






	

