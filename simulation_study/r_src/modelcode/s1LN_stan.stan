data {
	int<lower=0> n;  						// sample size
	int<lower=0> m; 						// number of sampled areas
	int<lower=0> M; 						// total number of areas
	vector[m] FPC;							// finite population correlation for sampling variance
	
	int<lower=0> q_u; 						// number of parameters
	matrix[n, q_u] x_u; 					// covariates
	
	int area_loc_id[m]; 					// id for when we change to a new area
	int area_loc_size[m]; 	    			// number of sampled individuals for each area
	
	int y[n];       						// outcome - int
	vector[n] y2;       					// outcome - real
	vector[n] w_tilde;						// full sample scaled weights
	int ps_area[n];							// area indicators
	vector[n] w;							// sample weights
	vector[n] w_ss2;						// area sample scaled weights squared
	vector[m] sum_w;						// sum of sample weights per area
	real bot_SR; 							// bottom of smoothing ratio (SR)
}
transformed data{
	// DECLARACTIONS
		vector[m] inv_sum_w2;					// sum of sample weights per area squared
		vector[n] w2; 							// sample weights squared
	
	// DEFINITIONS
		w2 = square(w);
		for(i in 1:m){
			inv_sum_w2[i] = pow(sum_w[i], -2.0);
		}
}
parameters {
	// LOGISTIC MODEL
	vector[q_u] beta_u;							// fixed parameters
	real intercept_u;							// intercept
	vector[m] Z_e;
	real<lower=0> sigma_e;
}
transformed parameters{
	vector[n] gmma;
	for(i in 1:n){
		gmma[i] = intercept_u + x_u[i,] * beta_u + Z_e[ps_area[i]] * sigma_e;
	}
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
		vector<lower=0>[m] gamma_pd; 					// sampling variance -> empirical logit
		vector<lower=0>[m] var_mu_pd; 					// sampling variance
		vector[M] mu_pd; 								// Pseudo-direct estimate 
		vector[m] theta_pd; 								// Pseudo-direct estimate -> empirical logit
		real log_lik[n];								// LogLikelihood 
		real SR;										// smoothing ratio (SR)
		
	// DEFINITIONS		
	{ // local scope 
			vector[n] p;								// p_ij
			vector[n] wp; 								// w_ij * p_ij		
			vector[n] pmy; 								// y_ij - p_ij							
			vector[n] wpmy; 							// w_ij * (y_ij - p_ij)
			vector[n] wy; 								// w_ij * y_ij
			vector[m] top_SR_per_area; 					// abs value of bias per area
			vector[m] bias_D;							// weighted mean of bias
			vector[m] mu_d;								// direct estimate
			
			// kernel of variance calculation
				vector[m] model_element;
				vector[m] data_element;
				vector[m] bias_element;
			p = inv_logit(gmma);
			
			// weights times probabilities
			wp = w .* p; 
			pmy = p - y2;
			wpmy = w .* pmy;
			wy = w .* y2;
			for(j in 1:m){
			
				// Weighted means
					mu_pd[j] = sum(segment(wp, area_loc_id[j], area_loc_size[j]))/sum_w[j];
					mu_d[j] = sum(segment(wy, area_loc_id[j], area_loc_size[j]))/sum_w[j];
					// sum( w * (p - y) ) / ( sum( w ) )
					bias_D[j] = sum(segment(wpmy, area_loc_id[j], area_loc_size[j]))/sum_w[j];
				
				// kernel of variance calculation
					// sum( w_ij^2 * (y - Direct)^2 )
					data_element[j] = sum( segment(w_ss2, area_loc_id[j], area_loc_size[j]) .* square(segment(y2, area_loc_id[j], area_loc_size[j]) - mu_d[j]) );
					// sum( w_ij^2 * (p - mu_pd)^2 )
					model_element[j] = sum( segment(w_ss2, area_loc_id[j], area_loc_size[j]) .* square(segment(p, area_loc_id[j], area_loc_size[j]) - mu_pd[j]) );
					// sum( w_ij^2 * ((p-y) - bias_D)^2 )
					bias_element[j] = sum( segment(w_ss2, area_loc_id[j], area_loc_size[j]) .* square(segment(pmy, area_loc_id[j], area_loc_size[j]) - bias_D[j]) );
				
				// Sampling variance for mu_pd
				var_mu_pd[j] = FPC[j] * ( data_element[j] + bias_element[j] );
				
				// Sampling variance for theta_pd
				theta_pd[j] = logit(mu_pd[j]);
				gamma_pd[j] = var_mu_pd[j] / ( mu_pd[j]^2 * ( 1 - mu_pd[j] )^2 );
				
				// get top of smoothing ratio (SR)
				top_SR_per_area[j] = fabs( bias_D[j] );
			}
			SR = sum( top_SR_per_area ) / bot_SR;
			
			// LogLikelihood
			for(i in 1:n)
				log_lik[i] = bernoulli_logit_lpmf( y[i] | gmma[i] );
	}  
}






