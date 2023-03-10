data {
	int<lower=0> m; 									// number of sampled areas
	int<lower=0> M; 									// total number of areas
	int<lower=0> m_s;									// number of stable areas: m_s < m
	int<lower=0, upper=1> anyunstable;
	int id_s[m_s];										// indexes for the stable areas
	int<lower=0> id_us[anyunstable ? m-m_s : 0];		// indexes for unstable 
	//int<lower=0> id_us[m-m_s];						// indexes for unstable 

	
	int<lower=0> q_a; 									// number of parameters
	matrix[M, q_a] z_a; 								// covariates 
	
	int<lower=0> q1_gvf;								// design matrix for gamma model (INCLUDING INTERCEPT)
	matrix[m, q1_gvf] T;									
	
	vector[m] theta_o; 									// observed theta
	vector<lower=0>[m_s] sqrt_gamma_o; 					// observed gamma
}
transformed data{
	// DECLARACTIONS
	matrix[M, q_a] Q_ast;
	matrix[q_a, q_a] R_ast;
	matrix[q_a, q_a] R_ast_inverse;
		
	// thin and scale the QR decomposition
	Q_ast = qr_thin_Q(z_a) * sqrt(M - 1);
	R_ast = qr_thin_R(z_a) / sqrt(M - 1);
	R_ast_inverse = inverse(R_ast);
}
parameters {
	// NORMAL MODEL
	vector[q_a] l_a;													// fixed parameters
	real intercept_a;													// intercept
	vector[M] Z_v;														// standard normal for error
	real<lower=0> sigma_v;												// standard deviation of error
	
	// gamma MODEL
	vector[q1_gvf] omega;												// fixed parameters
	real<lower=0> sigma_gvf;											// residual error
}
transformed parameters{
	vector[M] nu = intercept_a + Q_ast * l_a + Z_v * sigma_v;			// linear predictor
	vector[m] sqrt_gamma = rep_vector(1, m); 
	vector[m] sqrt_gamma_mu = T * omega; 								// sampling variances
	sqrt_gamma[id_s] = sqrt_gamma_o;
	sqrt_gamma[id_us] = exp( sqrt_gamma_mu[id_us] + 0.5 * square(sigma_gvf) );
}
model{
	// PRIORS
		// for NORMAL MODEL
			target += normal_lpdf( l_a | 0, 2 );
			target += student_t_lpdf(intercept_a | 3, 0, 1);
			target += std_normal_lpdf( Z_v );
			target += normal_lpdf( sigma_v | 0, 2 );
		// for gamma MODEL
			target += normal_lpdf( omega | 0, 2 );
			target += cauchy_lpdf( sigma_gvf | 0, 2);
		
	// GVFs
		target += normal_lpdf( log( sqrt_gamma[id_s] ) | sqrt_gamma_mu[id_s] , sigma_gvf );

	// MODEL	
		// NORMAL MODEL
		for(i in 1:m){
			target += normal_lpdf( theta_o[i] | nu[i], sqrt_gamma[i] );
		}
}
generated quantities{
	real log_lik[m];
	vector[q_a] lambda_a = R_ast_inverse * l_a;						// get correct regression coefficients
	// LARGE VALUES MEANS MORE OF THE DIRECT ESTIMATOR
	vector[m] comp_est_weight = sigma_v^2 ./ ( square(sqrt_gamma) + sigma_v^2 ); 	// composite estimator weight
	vector[M] mu = inv_logit(nu);
	for(i in 1:m)
		log_lik[i] = normal_lpdf( theta_o[i] | nu[i], sqrt_gamma[i] );
}






