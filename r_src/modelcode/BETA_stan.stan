functions {
	// function to transform sqrt_psi to unconstrained space
	vector f_sqrt_psi(vector x){
		return log( x ./ (0.5 - x));
	}
	// inverse function 
	vector f_sqrt_psi_inv(vector x){
		return ( 0.5 * exp(x) ) ./ ( 1 + exp(x) );
	}
}
data {
	int<lower=0> m; 									// number of sampled areas
	int<lower=0> M;										// total number of areas
	vector[M] logNi; 									// log population size in each area
	int<lower=0> m_s;									// number of stable areas: m_s < m
	int id_s[m_s];										// indexes for the stable areas
	int id_us_miss[M-m_s];								// indexes for unstable and missing  

	int<lower=0> q_a;           						// number of parameters
	matrix[M, q_a] z_a; 								// covariates
	
	vector<lower=0,upper=1>[m] mu_d_o; 					// observed mu_d_o
	vector<lower=0, upper=0.5>[m_s] sqrt_psi_o; 		// observed sqrt(psi)
	
	vector[2] logit_mu_bounds; 							// lower and upper bounds on logit mu
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
	// BETA MODEL
	real intercept_a;																// intercept
	vector[q_a] l_a;																// fixed parameters
	real<lower=0> sigma_v;															// variance of random effect
	vector[M] Z_v; 																	// standard normal variable
	
	// PSI MODEL
	real z_0;																		// fixed parameters
	real<upper=0> z_1;																// fixed parameters
	real<lower=0> sigma_sqrt_psi;														// residual error
}
transformed parameters{
// DECLARACTIONS
	vector[M] mu; 
	vector[M] eta;				
	vector[M] alpha;							
	vector[M] etaLower; 
	vector[M] etaUpper;
	vector[M] phi;
	vector[M] sqrt_psi = rep_vector(0, M); 
	vector[M] sqrt_psi_mu; 	
	vector[M] psi;
	
// create FULL VECTORS
	// GVF model
		sqrt_psi = rep_vector(0, M); 
		sqrt_psi_mu = z_0 + z_1 * logNi;
		sqrt_psi[id_s] = sqrt_psi_o;
		sqrt_psi[id_us_miss] = f_sqrt_psi_inv( sqrt_psi_mu[id_us_miss] );
		psi = square(sqrt_psi);
		
	// ensure bounds on mu are respected
	etaLower = logit(0.5 * (1 - sqrt(1 - 4 * psi)));
	etaUpper = logit(0.5 * (1 + sqrt(1 - 4 * psi)));
	for(i in 1:M){
		if (etaLower[i] < logit_mu_bounds[1])
			etaLower[i] = logit_mu_bounds[1];
		if (etaUpper[i] > logit_mu_bounds[2])
			etaUpper[i] = logit_mu_bounds[2];
	}
	// Non-centered parameterization: https://mc-stan.org/docs/2_28/stan-users-guide/reparameterization.html
	eta = intercept_a + Q_ast * l_a + Z_v * sigma_v;
	// adjust the bounds on eta and transform to mu: https://mc-stan.org/docs/2_28/stan-users-guide/vectors-with-varying-bounds.html
	mu = inv_logit(etaLower + (etaUpper - etaLower) .* inv_logit(eta));
	phi = ((mu .* (1 - mu)) ./ (psi)) - 1;
	alpha = mu .* phi;
	
}
model{
	// PRIORS
		// for BETA MODEL
			target += normal_lpdf(l_a | 0, 2);
			target += student_t_lpdf(intercept_a | 3, 0, 1);
			target += std_normal_lpdf( Z_v );
			target += std_normal_lpdf( sigma_v );
		// for PSI MODEL
			target += normal_lpdf( z_0 | 0, 2 );
			target += normal_lpdf( z_1 | 0, 2 );
			target += cauchy_lpdf( sigma_sqrt_psi | 0, 2 );
			
	// GVFs
		target += normal_lpdf( f_sqrt_psi( sqrt_psi[id_s] ) | sqrt_psi_mu[id_s] , sigma_sqrt_psi );

	// MODEL
		// Likelihood for stable mu_d_o's and phi's
		for(i in 1:m){
			target += beta_lpdf( mu_d_o[i] | alpha[i], phi[i] - alpha[i] );
		}
}
generated quantities{
	real log_lik[m];
	vector[q_a] lambda_a = R_ast_inverse * l_a;
	for(i in 1:m)
		log_lik[i] = beta_lpdf( mu_d_o[i] | alpha[i], phi[i] - alpha[i] );
}


