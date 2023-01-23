data {
	int<lower=0> its;								// number of draws from LOGISTIC MODEL
	int<lower=0> m; 								// number of sampled areas
	int<lower=0> m_s;								// number of stable areas: m_s < m
	int<lower=0> M; 								// total number of areas
	int<lower=0, upper=1> anyunstable;				// indicator for whether any areas are unstable
	int id_s[m_s];									// indexes for the stable areas
	int<lower=0> id_us[anyunstable ? m-m_s : 0];	// indexes for unstable 
	
	int<lower=0> q_a; 								// number of parameters
	matrix[M, q_a] z_a; 							// covariates 
	
	int<lower=0> q1_gvf;							
	matrix[m, q1_gvf] T;							// design matrix for gamma model (INCLUDING INTERCEPT)
	
	matrix[its, m] theta_pd_its;					// posterior draws of theta_pd from LOGISTIC MODEL
	vector<lower=0>[m] theta_pd_sd;					// posterior sd of theta_pd from LOGISTIC MODEL
	vector<lower=0>[m] sqrt_gamma_input;			// posterior mean of sqrt(gamma_pd) from LOGISTIC MODEL
}
transformed data{
	// DECLARACTIONS
	matrix[M, q_a] Q_ast;
	matrix[q_a, q_a] R_ast;
	matrix[q_a, q_a] R_ast_inverse;
	real<lower=0> PL_w = 1.0/its;
		
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
	vector[m] bar_theta_pd;												// column means of theta_pd_its
	
	// gamma MODEL
	vector[q1_gvf] omega;												// fixed parameters
	real<lower=0> sigma_gvf;											// residual error
}
transformed parameters{
	vector[M] nu = intercept_a + Q_ast * l_a + Z_v * sigma_v;
	vector[M] mu = inv_logit(nu);
	
	// GVF for sqrt_gamma
		vector[m] sqrt_gamma_mu = T * omega;
		vector[m] sqrt_gamma = sqrt_gamma_input;
		sqrt_gamma[id_us] = exp( sqrt_gamma_mu[id_us] + 0.5 * square(sigma_gvf ) );
}
model{
	// PRIORS
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
		for(i in 1:m){
			//target += normal_lpdf( theta_pd_its[,i] | bar_theta_pd[i], theta_pd_sd[i] );
			target += PL_w * normal_lpdf( theta_pd_its[,i] | bar_theta_pd[i], theta_pd_sd[i] );
			target += normal_lpdf( bar_theta_pd[i] | nu[i], sqrt_gamma[i] );
		}
}
generated quantities{
	// get correct regression coefficients
	vector[q_a] lambda_a = R_ast_inverse * l_a;	
	// composite estimator weight	
	vector[m] comp_est_weight = sigma_v^2 ./ ( square(sqrt_gamma) + sigma_v^2 ); 	// composite estimator weight
	vector[m] il_bar_theta_pd = inv_logit(bar_theta_pd);
}






