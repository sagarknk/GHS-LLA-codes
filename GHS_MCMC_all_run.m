n = 25;
q = 10;
burnin = 1e3;
nmc = 5e3;

prec_struc = 2; %%% can be 1 for 'hubs'

if prec_struc == 1
    struc_label = 'hubs';
else
    struc_label = 'random';
end

m = 5; %%% number of replicates 

for data_idx = 1:m

    fprintf("%d data is being processed\n", data_idx);
    %%%%%%%%%%
    X_mat = readmatrix(['./Data/GHS_sim_p',num2str(q),struc_label,num2str(n),'_data',num2str(data_idx),'.csv']);
    %%%%%%%%%%
    S = X_mat' * X_mat;
    %%%%%%%%%%
    GHS(S,n,burnin,nmc,data_idx,prec_struc)
end

