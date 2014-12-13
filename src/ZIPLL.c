
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/********************************************************************************************************
 invert inverts a function
 ********************************************************************************************************/

double invert(int n, double A[n][n], double Ainv[n][n])
{	
	// this is a copy of A to avoid overwrighting it.
	double X[n][n], Achol[n][n];
	int i, j, k;
	
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			Achol[i][j]=0;
			X[j][j]=0;
		}
	}
	
	//cholesky to get upper triangle matrix
	for(i=0; i<n; i++){
		Achol[i][i] = A[i][i];
		for(k=0; k<i; k++) Achol[i][i] -= Achol[k][i]*Achol[k][i];
		Achol[i][i] = sqrt(Achol[i][i]);
		
		for(j=i+1; j<n; j++){
			Achol[i][j] = A[i][j];
			for(k=0; k<i; k++) Achol[i][j] -= Achol[k][i]*Achol[k][j];
			Achol[i][j] /= Achol[i][i];
		}
		//also set to 0 this will be the inverse of the upper tri 
		for(j=0; j<n; j++){
			X[i][j] = 0;
		}
	}
	
	// invert from chol
	X[0][0] = 1/Achol[0][0];
	for(i=1; i<n; i++){   //column
		X[i][i] = 1/Achol[i][i];
		for(j=i-1; j>=0; j=j-1){  //row
			for(k=j+1; k<=i; k++) X[j][i] -= Achol[j][k]*X[k][i]/Achol[j][j];		
		}
	}
	
	
	//set Ainv 
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			Ainv[i][j] = 0;
			for(k=0; k<n; k++) Ainv[i][j] += X[i][k]*X[j][k];
		}
	}
	return 0;
}



 


/********************************************************************************************************
 cholesky decomposition.  not that t(Achol) %*% Achol = A, i.e. it is backwards from many decopositions.
 ********************************************************************************************************/



double chol(int n, double A[n][n], double Achol[n][n])
{
	int i,j,k;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			Achol[i][j]=0;
		}
	}
	for( i=0; i<n; i++){
		Achol[i][i] = A[i][i];
		for( k=0; k<i; k++) Achol[i][i] -= Achol[k][i]*Achol[k][i];
		Achol[i][i] = sqrt(Achol[i][i]);
		
		for( j=i+1; j<n; j++){
			Achol[i][j] = A[i][j];
			for( k=0; k<i; k++) Achol[i][j] -= Achol[k][i]*Achol[k][j];
			Achol[i][j] /= Achol[i][i];
		}
	}
	return 0;
}


/********************************************************************************************************
 generates a multivariate random normal domension n
 ********************************************************************************************************/

double rmvnorm(int n, double mu[], double Sigma[n][n], double x[])
{
	
	//generate standard normal
	double z[n];
	for(int i=0; i<n; i++){
		z[i] = rnorm(0,1);
	}
	
	//matrix for cholesky decomp to transform z
	double V[n][n];
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			V[i][j] = 0;
		}
	}
	chol(n, Sigma, V);
	
	
	//transform
	for(int i=0; i<n; i++){
		x[i] = mu[i];
		for(int j=0; j<n; j++){
			x[i] += z[j]*V[j][i];
		}
	}
	return 0;
}



/********************************************************************************************************
 generates an wishart random variable.  Note that this takes a variance matrix not a percision.
 ********************************************************************************************************/

void rwish_inv(int n , int nu, double Sigma[n][n], double W[n][n]){
	
	double  Z[n][n], ZC[n][n], CC[n][n], Sigma_inv[n][n];	
	int i, j, k, i2;

	// set everything to 0
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			CC[i][j]=0;
			ZC[i][j]=0;
			Z[i][j]=0;
			Sigma_inv[i][j]=0;
		}
	}
	//invert
	invert(n,Sigma,Sigma_inv);
    
	//cholesky
	for(i=0; i<n; i++){
		CC[i][i] = Sigma_inv[i][i];
		for(k=0; k<i; k++) CC[i][i] -= CC[k][i]*CC[k][i];
		CC[i][i] = sqrt(CC[i][i]);
		
		for(j=i+1; j<n; j++){
			CC[i][j] =Sigma_inv[i][j];
			for(k=0; k<i; k++) CC[i][j] -= CC[k][i]*CC[k][j];
			CC[i][j] /= CC[i][i];
		}
	}
	
	// start generating.
	for(i=0; i<n; i++){
		Z[i][i] = sqrt( rchisq(nu-i) ); 
		for(int j=i+1; j <n; j++){
			Z[i][j] = rnorm(0,1);
		}
	}
	
	//multiply the matricies
	for(i=0; i < n; i++){
		for(j=0; j < n; j++){
			ZC[i][j] = 0;
			for(i2=0; i2 < n; i2++){
				ZC[i][j] += Z[i][i2] * CC[i2][j];
			}
		}
	}
	
	//calculate XtX
	for(i=0; i < n; i++){
		for(j=0; j < n; j++){
			W[i][j] = 0;
			for(i2=0; i2 < n; i2++){
				W[i][j] += ZC[i2][i] * ZC[i2][j];
			}
		}
	}
	return;
}

		
				 
/********************************************************************************************************
rtnorm
This funtion draws a random variable from a truncated normal distribution
********************************************************************************************************/
				 
				 
				 
double rtnorm(double mu, double sigma,double trunpt, int above)
{
	/* function to draw truncated normal
	 above=1 means from above b=trunpt, a=-inf
	 above=0 means from below a=trunpt, b= +inf
	 */
	double tpt,rnd,result ;
				
	rnd=log(unif_rand());
	
	if (above) {		
		tpt=pnorm(((trunpt-mu)/(sigma)), 0.0, 1.0, 1, 1);
		result = mu + sigma*qnorm(rnd+tpt, 0.0, 1.0, 1, 1);
	}else {	
		tpt=pnorm(-((trunpt-mu)/(sigma)), 0.0, 1.0, 1, 1);
		result = mu - sigma*qnorm(rnd+tpt, 0.0, 1.0, 1, 1);
	}
				
	return result;
} 


/********************************************************************************************************
	g_func
	This funtion calculates G(x_{ijk},ac50_{ij},beta_{ij}) for ZIPLL.
********************************************************************************************************/

double g_func(double logdiff, double beta[], double knots[], int n_bases){
	double g=0;
	for(int l=1; l<=n_bases; l++){
		if(knots[l-1]<0 && logdiff<0){
			g += fmin( fmax(logdiff, knots[l-1]) - knots[l] ,0) * beta[l-1];
		}
		if(knots[l]>0 && logdiff>0){
			g += fmax( fmin(logdiff, knots[l]) - knots[l-1] ,0) * beta[l-1];
		}	
	}
	return g;
}



/********************************************************************************************************
 d_ss
 This computes the difference in sum of squares divided by 2 for two normal pdf. 
 ********************************************************************************************************/
double d_ss(int dim, double mu[], double Sigma[dim][dim], double y1[], double y2[]){
	
	double r1[dim], r2[dim], ss1=0, ss2=0;
	for(int i=0; i < dim; i++){
		r1[i] = y1[i] - mu[i];
		r2[i] = y2[i] - mu[i];
	}
	
	
	for(int j=0; j<dim; j++){
		for(int k=0; k<dim; k++){
			ss1 += r1[j]*r1[k]*Sigma[j][k];
			ss2 += r2[j]*r2[k]*Sigma[j][k];
		}
	}
	
	double diff = -(ss1-ss2)/2;
	
	return diff;
} 

double d_ss_vec(int dim, double mu[], double Sigma[], double y1[], double y2[])
{
	
	double r1[dim], r2[dim], ss1=0, ss2=0;
	for(int i=0; i < dim; i++){
		r1[i] = y1[i] - mu[i];
		r2[i] = y2[i] - mu[i];
	}
	
	
	for(int j=0; j<dim; j++){
		for(int k=0; k<dim; k++){
			ss1 += r1[j]*r1[k]*Sigma[j*dim+k];
			ss2 += r2[j]*r2[k]*Sigma[j*dim+k];
		}
	}
	
	double diff = -(ss1-ss2)/2;
	
	return diff;
} 


/********************************************************************************************************
 update_bottom
 This funtion updates the bottom
 ********************************************************************************************************/



void update_bottom(double y[], double x[], int n, double knots[], int n_bases,
				  double theta[], double theta_star[],
				  double beta_temp[], int z,
				  double Sigma_inv_mu[3][3], double mu[], double sig2_inv
				  ){
	
	
	
	double g = 0, g2 = 0, yg = 0, d_e_star_t = 0, d_y_b = 0, sumy = 0, 
		   E_GTT_b, V_GTT_b, E_LTT_b, V_LTT_b, GTT_log_prob_b, LTT_log_prob_b;
	int k;
	
	/*loop through observations and calculate needed values*/
	for(k=0; k<n; k++){
		sumy += y[k];
	}
	
	/* this is the top below the bottom case */
	V_GTT_b = 1/(n * sig2_inv + Sigma_inv_mu[1][1]);
	E_GTT_b = V_GTT_b * (  sig2_inv*sumy
						 - Sigma_inv_mu[0][1] * (theta_star[0] - mu[0])
						 + Sigma_inv_mu[1][1] *  mu[1]
						 - Sigma_inv_mu[2][1] * (theta_star[2] - mu[2]) );
	
	if(z==0){
		theta_star[1] = theta[1] = rnorm(E_GTT_b, sqrt(V_GTT_b));
	}else{
		
		for(k=0; k<n; k++){
			g = 1/(1+exp(g_func(log(x[k]) - log(theta[2]), beta_temp, knots, n_bases)));
			g2 += g * g;
			yg += (g) * (y[k] - theta[0] * (1-g));
			d_e_star_t += dnorm(y[k], theta[0]*(1-g), sqrt(1/sig2_inv), 1);
			d_y_b += dnorm(y[k],0,sqrt(1/sig2_inv),1);
		}
	
		/* for top greater than bottom case*/
		V_LTT_b = 1/(sig2_inv * g2 +  Sigma_inv_mu[1][1]);
		E_LTT_b = V_LTT_b * (  sig2_inv * yg
						 - Sigma_inv_mu[0][1] * (theta_star[0] - mu[0])
						 + Sigma_inv_mu[1][1] *  mu[1]
						 - Sigma_inv_mu[2][1] * (theta_star[2] - mu[2]) );
		
		/* calculate probablility below */
		LTT_log_prob_b = (pnorm(theta_star[0],E_LTT_b,sqrt(V_LTT_b), 1 , 1) 
						  - d_y_b
						  + d_e_star_t
						  - dnorm(0,E_LTT_b,sqrt(V_LTT_b), 1)
						  );
		GTT_log_prob_b = (pnorm(theta_star[0],E_GTT_b,sqrt(V_GTT_b), 0, 1) 
						  - dnorm(0,E_GTT_b,sqrt(V_GTT_b), 1) 
						  );
		double prob_GTT = 1 / (1+ exp(LTT_log_prob_b - GTT_log_prob_b));
		
		/*indicator for below*/
		
		double ru=unif_rand();
		if(ru>prob_GTT){
			theta_star[1] = theta[1] = rtnorm(E_LTT_b, sqrt(V_LTT_b), theta_star[0], 1);
		}else{
			theta_star[1] = theta[1] = rtnorm(E_GTT_b, sqrt(V_GTT_b), theta_star[0], 0);
		}
		
		
	}
} 


/********************************************************************************************************
 update_top
 This funtion updates the top
 ********************************************************************************************************/


void update_top(double y[], double x[], int n, double knots[], int n_bases,
				  double theta[], double theta_star[],
				  double beta_temp[], int z,
				  double Sigma_inv_mu[3][3], double mu[], double sig2_inv
				  ){

	
	double g = 0, g2 = 0, yg = 0, d_e_star_t = 0, d_y_b = 0, E_GTB_t, V_GTB_t,
		  E_LTB_t, V_LTB_t, GTB_log_prob_t, LTB_log_prob_t;
	int k;
	
	
	/* this is the top below the bottom case */
	V_LTB_t = 1/(Sigma_inv_mu[0][0]);
	E_LTB_t = V_LTB_t * (  Sigma_inv_mu[0][0] *  mu[0]
						 - Sigma_inv_mu[1][0] * (theta_star[1] - mu[1])
						 - Sigma_inv_mu[2][0] * (theta_star[2] - mu[2]) );
	
	if(z==0){
		theta_star[0] = rnorm(E_LTB_t, sqrt(V_LTB_t));
	}else{
	
		/*loop through observations and calculate needed values*/
		for(k=0; k<n; k++){
			g = 1/(1+exp(g_func(log(x[k]) - log(theta[2]), beta_temp, knots, n_bases)));
			g2 += (1-g) * (1-g);
			yg += (1-g) * (y[k] - theta[1]*g);
			d_e_star_t += dnorm(y[k], theta[1]*g, sqrt(1/sig2_inv), 1);
			d_y_b += dnorm(y[k],theta[1],sqrt(1/sig2_inv),1);
		}
		
		/* for top greater than bottom case*/
		V_GTB_t = 1/(sig2_inv * g2 +  Sigma_inv_mu[0][0]);
		E_GTB_t = V_GTB_t * (  sig2_inv * yg
							 + Sigma_inv_mu[0][0] *  mu[0]
						 - Sigma_inv_mu[1][0] * (theta_star[1] - mu[1])
						 - Sigma_inv_mu[2][0] * (theta_star[2] - mu[2]) );
						 
		/* calculate probablility below */
		GTB_log_prob_t = (pnorm(theta[1],E_GTB_t,sqrt(V_GTB_t), 0 , 1) 
						  - dnorm(0,E_GTB_t,sqrt(V_GTB_t), 1)
						  + d_e_star_t);
		LTB_log_prob_t = (pnorm(theta[1],E_LTB_t,sqrt(V_LTB_t), 1, 1) 
						  - dnorm(0,E_LTB_t,sqrt(V_LTB_t), 1) 
						  + d_y_b);
		double prob_below = 1 / (1+ exp(GTB_log_prob_t - LTB_log_prob_t));
	
		/*indicator for below*/
		double ru=unif_rand();
		if(ru>prob_below){
			theta_star[0] = rtnorm(E_GTB_t, sqrt(V_GTB_t), theta[1], 0);
		}else{
			theta_star[0] = rtnorm(E_LTB_t, sqrt(V_LTB_t), theta[1], 1);
		}
	}
	
	//update top
	theta[0] = fmax(theta_star[0],theta_star[1]);
} 




/********************************************************************************************************
 update_top
 This funtion updates the top
 ********************************************************************************************************/


int update_Z(double y[], double x[], int n, double knots[], int n_bases,
				  double theta[], int z,
				  double beta_temp[], 
				  double sig2_inv, double psi
				  ){
	
	double SS = 0, psi_hat;
	int k;
	
	for(k=0; k < n; k++){
		SS += (y[k] - theta[1])*(y[k] - theta[1])
			- (y[k] - theta[0] + (theta[0]-theta[1])/(1+exp(g_func(log(x[k]) - log(theta[2]), beta_temp, knots, n_bases))))
			* (y[k] - theta[0] + (theta[0]-theta[1])/(1+exp(g_func(log(x[k]) - log(theta[2]), beta_temp, knots, n_bases))));
	} 
	
	psi_hat = psi / (psi + (1-psi) * exp(-sig2_inv * SS / 2));
	
	double ru=unif_rand();
	if(ru<psi_hat){
		z = 1;
	}else{
		z = 0;
	}
	//Rprintf("ps_hat %f  %d\n", psi_hat, z);
	return z;
} 




/********************************************************************************************************
	tb_func
	This funtion updates theta, beta, and Z.
********************************************************************************************************/


void FIT_ZIPLL(double *resp, double *conc, int *dimx, int *n, 
				  double *knots, int *nknots,
				  double *BE_in, double *BEsd_in,
				  double *AC50_out, double *TOP_out, double *AC50sd_out, double *TOPsd_out,
				  double *active_out, double *kappa,
				  int *smax, int *burnin,
				  double *mu_0, double *Diag_Sigma_0
				  ){
	GetRNGstate();
	
	//-------------------------------------------------------------------------------------
	//indecies
	int i, i2, j, k, k2, u, chem, assay, rep;
	
	//-------------------------------------------------------------------------------------
	//find number of chemicals and assays
    int n_assays = dimx[1];
    int n_chems = dimx[0];
    int N = n_chems*n_assays;  //total number of curves
	int max_n = dimx[2];  //maximum number of obs allowed
	int n_star = 0;  //total number of observations across all curves
	for(i=0; i<N; i++) n_star += n[i];
	int K=*nknots-1;
	
		
	//-------------------------------------------------------------------------------------
	//Setup data storage

	
	double **x;
	x = (double **) calloc(N,sizeof(double*));
	for (i=0; i<N; i++) {
		x[i] = (double *) calloc(max_n, sizeof(double));
	}
	double **y;
	y = (double **) calloc(N,sizeof(double*));
	for (i=0; i<N; i++) {
		y[i] = (double *) calloc(max_n, sizeof(double));
	}	
	double **BE;
	BE = (double **) calloc(N,sizeof(double*));
	for (i=0; i<N; i++) {
		BE[i] = (double *) calloc(max_n, sizeof(double));
	}	
	double **BE2;
	BE2 = (double **) calloc(N,sizeof(double*));
	for (i=0; i<N; i++) {
		BE2[i] = (double *) calloc(max_n, sizeof(double));
	}
	
	
	
	for(i=0; i<N; i++){
		AC50_out[i] = 0;
		TOP_out[i] = 0;
		AC50sd_out[i] = 0;
		TOPsd_out[i] = 0;
		active_out[i] = 0;
		for(j=0; j<max_n; j++){
			x[i][j] = conc[j+i*dimx[2]];
			y[i][j] = resp[j+i*dimx[2]];
			BE[i][j] = 0; 
			BE2[i][j] = 0; 
		}
	}
	
	
	int Z[n_chems][n_assays];
	for(assay=0; assay<n_assays; assay++){
		for(chem=0; chem<n_chems; chem++){
			Z[chem][assay] = 1;
		}
    }

	//-------------------------------------------------------------------------------------
	//parameters, hyperparameters, and other prespecified numbers.
	
	double epsilon = 0.8, prblock = .75; //for resolvant
    double d1 = 1, d2 = 0.5; //for sig2
    double p1 = 1, p2 = 1; //for psi
    double g1 = 1, g2 = .5;  //for lambda
	
	//for mu and Sigma_mu
	double Sigma_0[3][3], Sigma_0_inv[3][3];
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			Sigma_0[i][j] = 0;
			Sigma_0_inv[i][j] = 0;
		}
	}
	Sigma_0[0][0] = Diag_Sigma_0[0];
	Sigma_0[1][1] = Diag_Sigma_0[1];
	Sigma_0[2][2] = Diag_Sigma_0[2];
	invert(3,Sigma_0,Sigma_0_inv);
	
    //for slope mean and variances
	int nu_mu = 6, nu_S=K+5;
    double V_mu[3][3], V_S[K][K];
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			V_mu[i][j] = 20*Sigma_0[i][j]*(nu_mu-5);
		}
	}
	for(k=0; k<K; k++){
		for(k2=0; k2<K; k2++){
			V_S[k][k2] = 0;
		}
		V_S[k][k] = 10;
	}
    double S0 = 2;

	
	//-------------------------------------------------------------------------------------
	// unknown parameters and assign starting values
		
	//Starting values
	double mu[3], // mean across assays
		   Sigma_bet_inv[3][3], // between assay inverse varinace
		   Sigma_S_inv[K][K], //inverse
		   S[K]
		   ;
		   
	//starting values for spline
	for(k=0; k<K; k++) S[k] = S0;
	for(i=0; i<3; i++) mu[i] = mu_0[i];
	
	
	
	double  t33[3][3],  m33[3][3], t3[3], m3[3], SSt_theta_a[3][3],
			sumbetastar[K], SSt_beta[K][K],
		    E_sk2[K], E_sk[K], V_sk_inv[K][K], V_sk[K][K];

		
	//Starting values for within assay variance
	for(assay=0; assay<n_assays; assay++){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				Sigma_bet_inv[i][j] = V_mu[i][j];
			}
			Sigma_bet_inv[i][i] = 1/V_mu[i][i];
		}
	}
	   
	rwish_inv(K, nu_S, V_S, Sigma_S_inv);
   

	double **beta_ac;
	beta_ac = (double **) calloc(N,sizeof(double*));
	for (i=0; i<N; i++) {
		beta_ac[i] = (double *) calloc(K, sizeof(double));
	}
	double **theta_ac;
	theta_ac = (double **) calloc(N,sizeof(double*));
	for (i=0; i<N; i++) {
		theta_ac[i] = (double *) calloc(3, sizeof(double));
	}
	double **theta_ac_star;
	theta_ac_star = (double **) calloc(N,sizeof(double*));
	for (i=0; i<N; i++) {
		theta_ac_star[i] = (double *) calloc(3, sizeof(double));
	}
	

		
	for(i=0; i<N; i++){
		rmvnorm(3, mu,Sigma_0,theta_ac_star[i]);
		for(k=0; k<K; k++){
			beta_ac[i][k] = exp(rnorm(S0,.1));
		}

		theta_ac[i][0] = fmax(theta_ac_star[i][0],theta_ac_star[i][1]);
		theta_ac[i][1] = theta_ac_star[i][1];
		theta_ac[i][2] = x[i][0] + (x[i][n[i]-1]-x[i][0]) / (1+exp(-theta_ac_star[i][2]));
	}
	

   
	double lambda = rgamma(d1, 1/d2);
	double sig2_inv = rgamma(g1, 1/g2);
	
	double psi = 0.5;
	
	
    //-------------------------------------------------------------------------------------
	//variables used throughout
	int n_resolv, row_temp, sumZ;
	double prop_AC50_star, prop_AC50, beta_star_temp[K], prop_beta_star[K], prop_beta[K], 
		   ran_ker, R, SSE = 0, y_hat, sumtheta_a[3], l_temp, theta_star_prop[3];
	

	//-------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------
	//MCMC	
	for(rep=0; rep<*smax; rep++){
	
		SSE = 0;  //sse for all curves
		
		//sum the assay specific parameters across all assays
		for(i=0; i<3; i++){
			sumtheta_a[i] = 0;
			for(i2=0; i2<3; i2++){
				SSt_theta_a[i2][i] = V_mu[i2][i] ;
			}
		}
	
		//sum the betas and squared betas across all chemicals and assays
		for(k=0; k<K; k++){
			sumbetastar[k] = 0;
			for(k2=0; k2<K; k2++){
				SSt_beta[k][k2] = V_S[k][k2] ;
			}
		}
		
		sumZ = 0;  //number of active curves in this assay
	
				
	//-------------------------------------------------------------------------------------								
	// assay loop	
		for(assay = 0; assay<n_assays; assay++){	
				
			
			
	//-------------------------------------------------------------------------------------
	// chemical loop
			for(chem = 0; chem < n_chems; chem++){
				row_temp = chem+assay*n_chems; // current row in the N length vector for this chem-assay combo
				
				
			
				//update bottom
				update_bottom(	y[row_temp], x[row_temp], n[row_temp], knots, K,
								theta_ac[row_temp], theta_ac_star[row_temp], 
								beta_ac[row_temp], Z[chem][assay],
								Sigma_bet_inv, mu, sig2_inv
								);
								
				//update top
				update_top(		y[row_temp], x[row_temp], n[row_temp], knots, K,
								theta_ac[row_temp], theta_ac_star[row_temp], 
								beta_ac[row_temp], Z[chem][assay],
								Sigma_bet_inv, mu, sig2_inv
								);
			
				//update Z
				 Z[chem][assay] = update_Z(	y[row_temp], x[row_temp], n[row_temp], knots, K,
								theta_ac[row_temp], Z[chem][assay],
								beta_ac[row_temp], sig2_inv, psi);
				
				
				/*number of updates*/
				n_resolv = trunc((log(unif_rand())/log(epsilon)));
			
				for(u=0; u<n_resolv; u++) {
					ran_ker = unif_rand();
								
					/* propose new value for AC50*/
					prop_AC50_star = theta_ac_star[row_temp][2] + rnorm(0, 1)*.5; 
					prop_AC50 = x[row_temp][0] + (x[row_temp][n[row_temp]-1]-x[row_temp][0]) / (1+exp(-prop_AC50_star));
				
					/* propose new values for Beta*/
					for(k=0; k<K; k++){
						beta_star_temp[k] = log(beta_ac[row_temp][k]);
						prop_beta_star[k] = beta_star_temp[k] + rnorm(0, 1)*.5;
						prop_beta[k] = exp(prop_beta_star[k]);
					}
				
					/* arrays of parms*/
					theta_star_prop[0] = theta_ac_star[row_temp][0];
					theta_star_prop[1] = theta_ac_star[row_temp][1];
					theta_star_prop[2] = prop_AC50_star;
				
				
					//Block update
					if(ran_ker<prblock){
						R = d_ss(3, mu, Sigma_bet_inv, theta_star_prop, theta_ac_star[row_temp]) + d_ss(K, S, Sigma_S_inv, prop_beta_star, beta_star_temp);
						if(Z[chem][assay]==1){
							for(k=0; k < n[row_temp]; k++){
								R += - sig2_inv * (
										  (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(prop_AC50), prop_beta, knots, K))))
										* (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(prop_AC50), prop_beta, knots, K)))) 
										- (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(theta_ac[row_temp][2]), beta_ac[row_temp], knots, K))))
										* (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(theta_ac[row_temp][2]), beta_ac[row_temp], knots, K))))
										)/2;
							}
						}
				
						if(R > log(unif_rand())){
							for(k=0; k<K; k++){ 
								beta_ac[row_temp][k] = prop_beta[k];
							}
							theta_ac_star[row_temp][2]= prop_AC50_star;
							theta_ac[row_temp][2] = prop_AC50;
						}
						
					}else{
					
						//theta first
						R = d_ss(3, mu, Sigma_bet_inv, theta_star_prop, theta_ac_star[row_temp]);
						if(Z[chem][assay]==1){
							for(k=0; k < n[row_temp]; k++){
								R += - sig2_inv * (
										  (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(prop_AC50), beta_ac[row_temp], knots, K))))
										* (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(prop_AC50), beta_ac[row_temp], knots, K)))) 
										- (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(theta_ac[row_temp][2]), beta_ac[row_temp], knots, K))))
										* (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(theta_ac[row_temp][2]), beta_ac[row_temp], knots, K))))
										)/2;
							}
						}
					
						if(R > log(unif_rand())){
							theta_ac_star[row_temp][2] = prop_AC50_star;
							theta_ac[row_temp][2] = prop_AC50;
						}
					
						//now beta
						// log ratio of multivariate normals
						R = d_ss(K, S, Sigma_S_inv, prop_beta_star, beta_star_temp);
						if(Z[chem][assay]==1){
							for(k=0; k < n[row_temp]; k++){
								R += - sig2_inv * (
										  (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(theta_ac[row_temp][2]), prop_beta, knots, K))))
										* (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(theta_ac[row_temp][2]), prop_beta, knots, K)))) 
										- (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(theta_ac[row_temp][2]), beta_ac[row_temp], knots, K))))
										* (y[row_temp][k] - theta_ac[row_temp][0] + (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][k]) - log(theta_ac[row_temp][2]), beta_ac[row_temp], knots, K))))
										)/2;		
							}
						}
					
					
						if(R > log(unif_rand())){
							for(k=0; k<K; k++){
								beta_ac[row_temp][k] = prop_beta[k]; 
							}
						}
					}
				}
			
			
			
		
					
				//update SSE and predicted values
				for(j=0; j<n[row_temp]; j++){
					// predicted value
					y_hat = (1-Z[chem][assay]) * (theta_ac[row_temp][1]) 
							 + Z[chem][assay]  * (theta_ac[row_temp][0] - (theta_ac[row_temp][0]-theta_ac[row_temp][1])/(1+exp(g_func(log(x[row_temp][j]) - log(theta_ac[row_temp][2]), beta_ac[row_temp], knots, K))));
					//predict at current locations
					if(rep>*burnin-1){
						BE[row_temp][j] += y_hat/(*smax-*burnin);  
						BE2[row_temp][j] += y_hat*y_hat/(*smax-*burnin);  
					}
					//records SSE
					SSE += ( y_hat - y[row_temp][j]) * ( y_hat - y[row_temp][j]);  // SSE
				}
				
				
				
				//keep sum for updating theta_a
				sumtheta_a[0] += theta_ac_star[row_temp][0];
				sumtheta_a[1] += theta_ac_star[row_temp][1];
				sumtheta_a[2] += theta_ac_star[row_temp][2];	   
				//sum of square parameters for within assay variance update
				SSt_theta_a[0][0] += (theta_ac_star[row_temp][0] - mu[0]) * (theta_ac_star[row_temp][0] - mu[0]);
				SSt_theta_a[0][1] += (theta_ac_star[row_temp][0] - mu[0]) * (theta_ac_star[row_temp][1] - mu[1]);
				SSt_theta_a[0][2] += (theta_ac_star[row_temp][0] - mu[0]) * (theta_ac_star[row_temp][2] - mu[2]);
				SSt_theta_a[1][0] += (theta_ac_star[row_temp][1] - mu[1]) * (theta_ac_star[row_temp][0] - mu[0]);
				SSt_theta_a[1][1] += (theta_ac_star[row_temp][1] - mu[1]) * (theta_ac_star[row_temp][1] - mu[1]);
				SSt_theta_a[1][2] += (theta_ac_star[row_temp][1] - mu[1]) * (theta_ac_star[row_temp][2] - mu[2]);
				SSt_theta_a[2][0] += (theta_ac_star[row_temp][2] - mu[2]) * (theta_ac_star[row_temp][0] - mu[0]);
				SSt_theta_a[2][1] += (theta_ac_star[row_temp][2] - mu[2]) * (theta_ac_star[row_temp][1] - mu[1]);
				SSt_theta_a[2][2] += (theta_ac_star[row_temp][2] - mu[2]) * (theta_ac_star[row_temp][2] - mu[2]);
			
				//sum Z
				sumZ += Z[chem][assay];
							
				//sum of squares for S
				for(k=0; k<K; k++){
						sumbetastar[k] += log(beta_ac[row_temp][k]);
					for(k2=0; k2<K; k2++){
						SSt_beta[k][k2] += (log(beta_ac[row_temp][k])-S[k]) * (log(beta_ac[row_temp][k2])-S[k2]) ;
					}
				}
			
						
			if(rep>*burnin-1){
				TOP_out[row_temp] += theta_ac[row_temp][0]/(*smax-*burnin);
				AC50_out[row_temp] += theta_ac[row_temp][2]/(*smax-*burnin);
				TOPsd_out[row_temp] += theta_ac[row_temp][0]*theta_ac[row_temp][0]/(*smax-*burnin);
				AC50sd_out[row_temp] += theta_ac[row_temp][2]*theta_ac[row_temp][2]/(*smax-*burnin);
				if(Z[chem][assay]==1 && (theta_ac[row_temp][0]-theta_ac[row_temp][1]) > *kappa) active_out[row_temp] += (double) 1/(*smax-*burnin);
			}
																
			}// end chemical loop
			
			//update psi for each assay
			psi = rbeta(p1 + sumZ, N+p2-sumZ);
			
			
		}// end assay loop
		//-------------------------------------------------------------------------------------
	
		
		
		//update the variance
		sig2_inv = rgamma((d1 + n_star/2), (1/(d2 +  SSE/2)));
		
		//update lambda
		l_temp = (S0 - S[0])*(S0 - S[0]);
		for(k=1; k<K-1; k++) l_temp += (S[k] - S[k-1])*(S[k] - S[k-1]);
		lambda = rgamma(g1 + K/2, 1/( g2 + l_temp/2));
		
		//update the between assay variance		
		rwish_inv(3, N+nu_mu, SSt_theta_a, Sigma_bet_inv);
		
		
		//update the overall mu
		for(i=0; i<3; i++){
			t3[i] = 0;
			for(i2=0; i2<3; i2++){
				t33[i][i2] = N*Sigma_bet_inv[i][i2] + Sigma_0_inv[i][i2];
				t3[i] += Sigma_bet_inv[i][i2] * sumtheta_a[i2] + Sigma_0_inv[i][i2] * mu_0[i2];
			}
		}
		
		invert(3,t33,m33);
		
		for(i=0; i<3; i++){
			m3[i] = 0;
			for(i2= 0; i2<3; i2++){
				m3[i] += m33[i][i2] * t3[i2];   
			}				   
		}
		
		rmvnorm(3, m3, m33, mu);
	

		//update the variance about the S spline
		rwish_inv(K, N+nu_S, SSt_beta, Sigma_S_inv);
		
		
		//variances for each S
		for(k=0; k<K; k++){
			//variance
			for(k2=0; k2<K; k2++){
				V_sk_inv[k][k2] = N * Sigma_S_inv[k][k2];
			}
			if(k<K-1){
				V_sk_inv[k][k] += 2*lambda;
			}else{
				V_sk_inv[k][k] += lambda;
			}
			//mean
			E_sk[k] = 0;
			if(k>0) E_sk[k] += S[k-1];
			if(k<K-1) E_sk[k] += S[k+1]; 
			if(k==0) E_sk[k] += S0;
			E_sk[k] = E_sk[k]*lambda;
		}
		invert(K, V_sk_inv,V_sk);
		
		for(k=0; k<K; k++){
			for(k2=0; k2<K; k2++){
				E_sk[k] += Sigma_S_inv[k][k2]*sumbetastar[k2];
			}
		}
		for(k=0; k<K; k++){
			E_sk2[k] = 0;
			for(k2=0; k2<K; k2++){
			E_sk2[k] += V_sk[k][k2]*E_sk[k2];
			}
		}

	
		//update matrix of S
		rmvnorm(K, E_sk2, V_sk, S);

		
	
	
	} // end MCMC
	//-------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------				
	
	
		
	//-------------------------------------------------------------------------------------
	for(i=0; i<N; i++){
		TOPsd_out[i] = sqrt(TOPsd_out[i] - TOP_out[i]*TOP_out[i]); 
		AC50sd_out[i] = sqrt(AC50sd_out[i] - AC50_out[i]*AC50_out[i]); 
		for(j=0; j<max_n; j++){
			BE_in[j+i*max_n] = BE[i][j];
			BEsd_in[j+i*max_n] = sqrt(BE2[i][j] - BE[i][j]*BE[i][j]);
		}
	}
	
		
	
	//-------------------------------------------------------------------------------------				
	// clean memory
	
	for(i = 0; i < N; i++){  
       free(beta_ac[i]);  
       free(theta_ac[i]);  
	   free(theta_ac_star[i]); 
       free(x[i]);  
       free(y[i]);  
       free(BE[i]); 
       free(BE2[i]);    
    }  
    free(beta_ac); 
    free(theta_ac); 
    free(theta_ac_star);
    free(x); 
    free(y); 
    free(BE); 
    free(BE2); 
	
		
	PutRNGstate();
	return;
}
