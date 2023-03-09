% Yunfan Li, December 2017
% Sampling for graphical horseshoe

function  GHS(S,n,burnin,nmc,data_idx,prec_struc)
% GHS MCMC sampler using data-augmented
% block (column-wise) Gibbs sampler
% Input:
%     S = Y'*Y : sample covariance matrix * n
%     n: sample size
%     burnin, nmc : number of MCMC burnins and saved samples

if prec_struc == 1
    struc_label = 'hubs';
else
    struc_label = 'random';
end
%%%%%%%%%%

[p] = size(S,1); indmx = reshape([1:p^2],p,p); 
upperind = indmx(triu(indmx,1)>0); 

indmx_t = indmx';
lowerind = indmx_t(triu(indmx_t,1)>0); 

omega_save = zeros(p,p,nmc);
lambda_sq_save = zeros(p*(p-1)/2,nmc);
omega_vector_save = zeros(0.5*p*(p-1), nmc);
tau_sq_save = zeros(1, nmc);

ind_noi_all = zeros(p-1,p);
for i = 1:p
       if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       
       ind_noi_all(:,i) = ind_noi;
end

% set initial values
Omega = eye(p); Sigma = eye(p);
Lambda_sq(1:p,1:p) = 1; Nu(1:p,1:p) = 1; tau_sq = 1; xi = 1;

tic;
for iter = 1: burnin+nmc    
          
   if(mod(iter,1000)==0)
    fprintf('iter = %d \n',iter);
   end
       
%%% sample Sigma and Omega=inv(Sigma)
    for i = 1:p
      ind_noi = ind_noi_all(:,i);     
      Sigma_11 = Sigma(ind_noi,ind_noi); sigma_12 = Sigma(ind_noi,i);
      sigma_22 = Sigma(i,i);
      s_21 = S(ind_noi,i); s_22 = S(i,i);
      lambda_sq_12 = Lambda_sq(ind_noi,i); nu_12 = Nu(ind_noi,i);
      %% sample gamma and beta
      gamma = gamrnd((n/2+1),2/s_22);    % random gamma with shape=n/2+1, rate=s_22/2
      inv_Omega_11 = Sigma_11 - sigma_12*sigma_12'/sigma_22;
      inv_C = s_22*inv_Omega_11+diag(1./(lambda_sq_12*tau_sq));
      inv_C_chol = chol(inv_C);
      mu_i = -inv_C\s_21;
      beta = mu_i+ inv_C_chol\randn(p-1,1);
      omega_12 = beta; omega_22 = gamma + beta'*inv_Omega_11*beta;
      %% sample lambda_sq and nu
      rate = omega_12.^2/(2*tau_sq)+1./nu_12;
      lambda_sq_12 = 1./gamrnd(1,1./rate);    % random inv gamma with shape=1, rate=rate
      nu_12 = 1./gamrnd(1,1./(1+1./lambda_sq_12));    % random inv gamma with shape=1, rate=1+1/lambda_sq_12
      %% update Omega, Sigma, Lambda_sq, Nu
      Omega(i,ind_noi) = omega_12; Omega(ind_noi,i) = omega_12;
      Omega(i,i) = omega_22;
      temp = inv_Omega_11*beta;
      Sigma_11 = inv_Omega_11 + temp*temp'/gamma;
      sigma_12 = -temp/gamma; sigma_22 = 1/gamma;
      Sigma(ind_noi,ind_noi) = Sigma_11; Sigma(i,i) = sigma_22;
      Sigma(i,ind_noi) = sigma_12; Sigma(ind_noi,i) = sigma_12;
      Lambda_sq(i,ind_noi) = lambda_sq_12; Lambda_sq(ind_noi,i) = lambda_sq_12;
      Nu(i,ind_noi) = nu_12; Nu(ind_noi,i) = nu_12;
    end
    
%%% sample tau_sq and xi
    omega_vector = Omega(tril(true(size(Omega)),-1));
    lambda_sq_vector = Lambda_sq(tril(true(size(Lambda_sq)),-1));
    rate = 1/xi + sum(omega_vector.^2./(2*lambda_sq_vector));
    tau_sq = 1/gamrnd((p*(p-1)/2+1)/2, 1/rate);    % inv gamma w/ shape=(p*(p-1)/2+1)/2, rate=rate
    xi = 1/gamrnd(1,1/(1+1/tau_sq));    % inv gamma w/ shape=1, rate=1+1/tau_sq

%%% save Omega, lambda_sq, tau_sq
    if iter >burnin           
         omega_save(:,:,iter-burnin) = Omega;
		 omega_vector_save(:,iter-burnin) = omega_vector;
         lambda_sq_save(:,iter-burnin) = lambda_sq_vector;
         tau_sq_save(1, iter-burnin) = tau_sq;
    end

end
time_taken = toc;
omega_save_final = mean(omega_save,3);

FileName=['./Results/GHS_MCMC_Workspace_of_',num2str(data_idx),'st_data_set_',num2str(n),'_',num2str(p),struc_label,'.mat'];

save(FileName, 'omega_save_final','omega_vector_save','time_taken','tau_sq_save','lambda_sq_save');
end

