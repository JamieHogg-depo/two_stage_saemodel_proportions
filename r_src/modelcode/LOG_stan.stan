data {
	int<lower=0> n;  						// sample size
	int<lower=0> m; 						// number of sampled areas
	int<lower=0> M; 						// total number of areas
	int<lower=0> N; 						// population size
	
	int<lower=0> q_u; 						// number of parameters
	matrix[n, q_u] x_u; 					// covariates
	matrix[N, q_u] x_u_census;				// covariates for census
	
	int area_loc_id[M]; 					// id for when we change to a new area
	int area_loc_size[M]; 	    			// number of sampled individuals for each area
	
	int y[n];       						// outcome
	vector[n] w_tilde;
	int ps_area[n];							// area indicators
	int ps_area_census[N];					// area indicators for census
}
parameters {
	// LOGISTIC MODEL
	vector[q_u] beta_u;							// fixed parameters
	real intercept_u;							// intercept
	vector[M] Z_e;
	real<lower=0> sigma_e;
}
transformed parameters{
	vector[n] gmma;
	for(i in 1:n)
		gmma[i] = intercept_u + x_u[i,] * beta_u + Z_e[ps_area[i]] * sigma_e;
}
model{
	// PRIORS
		// LOGISTIC MODEL
		target += normal_lpdf( beta_u | 0, 2 );
		target += std_normal_lpdf( Z_e );
		target += std_normal_lpdf( sigma_e );
		target += student_t_lpdf( intercept_u | 3, 0, 1 );
			
	// LOGISTIC MODEL
	for(i in 1:n)
		target += w_tilde[i] * bernoulli_logit_lpmf( y[i] | gmma[i] );
}
generated quantities{
	// DECLARACTIONS
		vector[M] mu; 		// area proportions						
		real log_lik[n];	// LogLikelihood 
		
	// DEFINITIONS		
	{ // local scope 
			vector[N] il_gmma_census;
			vector[N] Z_e_full;
			for(k in 1:N)
				Z_e_full[k] = Z_e[ps_area_census[k]];
			il_gmma_census = inv_logit( intercept_u + x_u_census * beta_u + Z_e_full * sigma_e );
			for(j in 1:M)
				mu[j] = mean(segment(il_gmma_census, area_loc_id[j], area_loc_size[j]));
			for(i in 1:n)
				log_lik[i] = bernoulli_logit_lpmf( y[i] | gmma[i] );
	}  
}






