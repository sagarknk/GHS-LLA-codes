total_sets = 1;

%n_EMs = 1;
%n_EMs = 10;
%n_EMs = 20;
n_EMs = 50;
%n_EMs = 100;

% starting point generation

Omega_saves = zeros(100,100,n_EMs);
rng(123456789);

for i = 1:n_EMs

start_point = eye(100);

    for row = 2:100
        d = 0;
            while d~=100
                row_seq = row:1:100;
                col_seq = 1:1:(100-row+1);

                rand_noise = -0.05 + rand(1, length(row_seq))*2*0.05;
                %rand_noise = -0.1 + rand(1, length(row_seq))*2*0.1;
                
                lin_idcs = sub2ind(size(start_point), row_seq, col_seq);
                start_point(lin_idcs) = rand_noise;

                lin_idcs = sub2ind(size(start_point), col_seq, row_seq);
                start_point(lin_idcs) = rand_noise;

                d = eig(start_point);
                d = sum(d>0);

            end 
    end 
    
fprintf("Finished %d data set generation out of %d data sets \n", i, n_EMs);
Omega_saves(:,:,i) = start_point;
end 

matObj = matfile("HS_LLA_LAPLACE_mix.mat");
dawson_vals = matObj.dawson_vals;
U_grid_linear = matObj.U_grid_linear;

%%%%%%%%%%%
dawson_vals = dawson_vals(1:3:length(dawson_vals));
U_grid_linear = U_grid_linear(1:3:length(U_grid_linear));
%%%%%%%%%%%
%fine_lambda_grid = [1e-6:1e-6:5e-1];
fine_lambda_grid = [1e-6:1e-5:5e-1];
%coarse_lambda_grid = [0.51:0.01:2000];
coarse_lambda_grid = [0.5:0.1:2000];

fun_cauchy_num = @(lambda) (1./(lambda.^3)).*(1./(1+lambda.*lambda));
fun_cauchy_denom = @(lambda) (1./(lambda)).*(1./(1+lambda.*lambda));

fine_cauchy_num_vals = fun_cauchy_num(fine_lambda_grid);
fine_cauchy_denom_vals = fun_cauchy_denom(fine_lambda_grid);

coarse_cauchy_num_vals = fun_cauchy_num(coarse_lambda_grid);
coarse_cauchy_denom_vals = fun_cauchy_denom(coarse_lambda_grid);

%%%%%%%%%%%%%
tau_val = 1/sqrt(120);
%%%%%%%%%%%%%
    
fine_lambda_tau_sq_inv = (1/tau_val^2).*(1./fine_lambda_grid)...
                        .*(1./fine_lambda_grid);

coarse_lambda_tau_sq_inv = (1/tau_val^2).*(1./coarse_lambda_grid)...
                        .*(1./coarse_lambda_grid);
      
for file_iter = 1:total_sets
    
    fprintf('Data set number - %d is being processed \n ',file_iter);
    FileName=['~/Documents/New_EM_n120_p100/New_Hubs_p100_n120/GHS_sim_p100hubs120_data',num2str(file_iter),'.csv'];
    xx = csvread(FileName);
    [n,q] = size(xx);
    S = xx'*xx;
    
    [Omega_est,total_iterations, each_time_taken] = ...
        Multi_start_point_Fixed_tau_HS_LLA_Laplace(Omega_saves, S, n, q, n_EMs,...
        U_grid_linear,dawson_vals,tau_val);
    
    FileName=['~/Documents/New_EM_n120_p100/New_Hubs_p100_n120/HS_LLA_MATs/HS_LLA_laplace_mix_Workspace_of_',num2str(file_iter),'st_data_set_with_',num2str(n_EMs),...
        '_start_points','.mat'];
    
    save(FileName, 'Omega_est', 'total_iterations',...
        'each_time_taken');
    %%%%%%%%%%%%%
    [Omega_est,total_iterations, each_time_taken] = ...
        Multi_start_point_Fixed_tau_HS_LLA_Cauchy(Omega_saves, S, n, q, n_EMs,...
        fine_cauchy_num_vals, coarse_cauchy_num_vals,...
        fine_cauchy_denom_vals, coarse_cauchy_denom_vals,...
        fine_lambda_tau_sq_inv, coarse_lambda_tau_sq_inv, tau_val);
    
    FileName=['~/Documents/New_EM_n120_p100/New_Hubs_p100_n120/HS_LLA_MATs/HS_LLA_cauchy_mix_Workspace_of_',num2str(file_iter),'st_data_set_with_',num2str(n_EMs),...
            '_start_points','.mat'];

    save(FileName, 'Omega_est', 'total_iterations',...
            'each_time_taken'); 
    %%%%%%%%%%%%%
    fprintf('Data set number - %d is finished \n ',file_iter);
    
end

% fun_exp_mixture= @(u,x) u.*sign(x).*exp(-abs(x).*u).*dawson(u./sqrt(2));
% fun_exp_mixture_derivative=@(u,x) exp(-abs(x).*u).*dawson(u./sqrt(2));
% 
% x = [-1:0.01:1];
% y = [0:0.5:100];
% 
% integral_1 = zeros(201,201);
% integral_2 = zeros(201,201);
% 
% for i = 1:201
%     i
%     for j = 1:201
%         if mod(j, 50) == 0
%             fprintf("%d J done\n", j)
%         end
%         integral_1(i,j) =  integral(@(u) fun_exp_mixture_derivative(u,x(i)),0,y(j));
%         integral_2(i,j) =  integral(@(u) fun_exp_mixture(u,x(i)),0,y(j));
%     end
% end 
% 
% G_B_num = integral(@(u) fun_exp_mixture_derivative(u,B),0,100);
% G_B_denom = integral(@(u)fun_exp_mixture(u,B),0,100);
%          
% copy_integral_1 = integral_1;
% copy_integral_2 = integral_2;
% copy_integral_1 = round(copy_integral_1,4);
% copy_integral_2 = round(copy_integral_2,4);
% 
% diff_copy_integral_1 = copy_integral_1(:,2:201) - copy_integral_1(:,1:200);
% diff_copy_integral_2 = copy_integral_2(:,2:201) - copy_integral_2(:,1:200);
% 
% diff_sum_1 = sum(diff_copy_integral_1.*diff_copy_integral_1,1);
% diff_sum_2 = sum(diff_copy_integral_2.*diff_copy_integral_2,1);
% 
% min_index_1 = find(diff_sum_1==min(diff_sum_1));
% min_index_2 = find(diff_sum_2==min(diff_sum_2));
% 
% 
%     






