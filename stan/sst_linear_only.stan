functions{
  //function to build the range matrix
  matrix inv_phi_mat_sq(real Phi_x, real Phi_y){
    vector[2] phi; 
    matrix[2,2] D;
      
      phi[1] =  Phi_x;
      phi[2] =  Phi_y;
      
      D[1,1] = 1/(2*phi[1]^2);
      D[1,2] = 0;
      D[2,1] = 0;
      D[2,2] = 1/(2*phi[2]^2);
    return D;
  }

  //function to estimate the distance covariance matrix
  matrix build_corr(real  Phi_x, real Phi_y, array[,] vector Dist,int N){
    matrix[2,2] Range;
    matrix[N,N] K;
    matrix[N,N] L;
    
    Range = inv_phi_mat_sq(Phi_x, Phi_y);

    for (i in 1:(N-1)){
      for (j in (i+1):N) {
        K[i,j] = exp(-quad_form(Range,Dist[i,j]));
        K[j,i] = K[i,j];
    }}
    for (k in 1:N){K[k,k] = 1;}
    
    L = cholesky_decompose(K);  
    return L;
  }

//function to estimate the kronecker product of two matrix
matrix kronecker_prod(matrix A, matrix B) {
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
}

//function to creata a matrix from the K product
matrix kronecker_output(matrix A, matrix B,vector E, int N,int P){
  matrix[rows(A) * rows(B), cols(A) * cols(B)] L_Kronecker;
  matrix[N,P] out; 

  L_Kronecker = kronecker_prod(A,B);
  out = to_matrix(L_Kronecker *E, N, P, 0);
  
  return(out);
}

//function to get the sd of params
vector extract_diag(matrix A){
  matrix[rows(A),cols(A)] B;
  vector[rows(A)] sigma; 
  
  B = multiply_lower_tri_self_transpose(A);
  for (i in 1:rows(A)){
      sigma[i] = sqrt(B[i,i]);
  }
  
  return(sigma);
}

}
data{
  //Indexes
  int<lower=1> NYear;
  int<lower=1> MAxYear;
  int<lower=1> NRegions;
  int<lower=1> NParams;
  //Observations
  matrix[NYear,NRegions] Y;
  matrix[NYear,NRegions] sigma_obs;
  matrix[NYear,NRegions] var_obs;
  vector[NYear] Trend;
  //Distance matrix
  array[NRegions, NRegions] vector[2] distance_vec; 
}
parameters {
  //Population parameters
  vector[NParams] mu_param;
  vector<lower=0>[NRegions] sigma;
  //SVP hyper parameters
  real<lower=0> phi_x;
  real<lower=0> phi_y;
  cholesky_factor_corr[NParams] L_rho;
  vector<lower=0>[NParams] L_xi;
  vector[NParams*NRegions] param_epsilon_raw;
  //Growth error
  matrix[MAxYear,NRegions] epsilon_N_raw;
}
transformed parameters{
    //SVP matrix
  matrix[NRegions,NRegions] L_Kappa;                        
  matrix[NParams,NParams] L_Tau;                                                                          
  matrix[NRegions,NParams] param_epsilon;
  //SVP parameters
  vector[NRegions] alpha;
  vector[NRegions] beta1;
  //Population
  matrix[MAxYear,NRegions] epsilon_N;
  matrix[MAxYear,NRegions] log_N;
  
  //spatial covariance
  L_Kappa =  build_corr(phi_x,phi_y,distance_vec,NRegions); 
  L_Tau = diag_pre_multiply(L_xi,L_rho);
  //Kronecker product
  param_epsilon = kronecker_output(L_Kappa,L_Tau,param_epsilon_raw,NRegions,NParams);
  
  //slice parameters 
  alpha = mu_param[1] + to_vector(param_epsilon[,1]);
  beta1 = mu_param[2] + to_vector(param_epsilon[,2]);
  
  //estimate the growth error 
  for(n in 1:NRegions){
    for(t in 1:MAxYear){
      epsilon_N[t,n] = sigma[n] * epsilon_N_raw[t,n];
  }}  
  
  
  //Population process
  for(n in 1:NRegions){
    for(t in 1:MAxYear){
      log_N[t,n] = alpha[n] + beta1[n]*Trend[t] + epsilon_N[t,n];
  }}
}
model{
  //Population parameter (alpha,beta)
  target += normal_lpdf(mu_param|0,0.25);
  target += gamma_lpdf(phi_x|5,0.25);
  target += gamma_lpdf(phi_y|5,0.25);
  target += lkj_corr_cholesky_lpdf(L_rho|2);
  target += student_t_lpdf(L_xi|4,0,0.5);
  target += normal_lpdf(param_epsilon_raw|0,1);
  
  //process error 
  target += student_t_lpdf(sigma|4,0,1);
  for(n in 1:NRegions){
    for (t in 1:MAxYear){
      target += normal_lpdf(epsilon_N_raw[t,n]|0,1);
  }}

  //Observation process
  for(n in 1:NRegions){
    for(t in 1:MAxYear){
      target += lognormal_lpdf(Y[t,n] |log_N[t,n], sigma_obs[t,n]);
    }}
}
generated quantities{
  //derived parameters
  matrix[NParams,NParams] rho;
  vector[NParams] sigma_param;
  //Loglikelihhod 
  matrix[MAxYear,NRegions] log_lik;
  //Prediction at t+1
  vector[NRegions] epsilon_pred_raw;
  vector[NRegions] epsilon_pred;
  vector[NRegions] pred_log_N;
  vector[NRegions] lik_pred; 

  //derived parameters
  rho = multiply_lower_tri_self_transpose(L_rho);
  sigma_param = extract_diag(L_Tau);
  
  
  //log likelihood 
  for(n in 1:NRegions){
    for(t in 1:MAxYear){
      log_lik[t,n] =  lognormal_lpdf(Y[t,n] |log_N[t,n], sigma_obs[t,n]);
  }}
  
  //Forecast next year and test fit
  for(n in 1:NRegions){
    epsilon_pred_raw[n] = normal_rng(0, 1);
    epsilon_pred[n] = sigma[n] * epsilon_pred_raw[n];
  }
  
  for(n in 1:NRegions){
     pred_log_N[n] =  alpha[n] + beta1[n]*Trend[(MAxYear+1)] + epsilon_pred[n];
  }
  
  for(n in 1:NRegions){
       lik_pred[n] =  lognormal_lpdf(Y[(MAxYear+1),n] |pred_log_N[n], sigma_obs[(MAxYear+1),n]);
  }
}
