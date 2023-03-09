n = 25;
q = 10;

n_EMs = 50; %%% number of start points for computing avg. results for LLA (l) or LLA (c)

prec_struc = 2; %%% can be 1 for 'hubs'

m = 5; %%% number of replicates 

for data_idx = 1:m
    CV_HS_LLA_Laplace(data_idx, n, q, prec_struc)
    CV_HS_LLA_Cauchy(data_idx, n, q, prec_struc)
end

for data_idx = 1:m
    Multi_start_point_Fixed_tau_HS_LLA_Laplace(data_idx, n_EMs, n, q, prec_struc)
    Multi_start_point_Fixed_tau_HS_LLA_Cauchy(data_idx, n_EMs, n, q, prec_struc)
end
    