
file_index = 1:5; %%% number of replications 

n = 25;
q = 10;

n_EMs = 50; %%% number of start points for computing avg. results for LLA (l) or LLA (c)

mix_type = 'laplace' ; %%% can be 'cauchy'
prec_struc = 2; %%% can be 1 for 'hubs'
%%%%%%%%%%
if prec_struc == 1
    struc_label = 'hubs';
else
    struc_label = 'random';
end
%%%%%%%%%%


Omega_each_DS = zeros(q,q,length(file_index));

for i = 1:length(file_index)

    file_iter = file_index(i);
    FileName=['./Results/HS_LLA_',mix_type,'_mix_Workspace_of_',num2str(file_iter),'st_data_set_with_',num2str(n_EMs),...
        '_start_points_',num2str(n),'_',num2str(q),struc_label,'.mat'];

    matObj = matfile(FileName);

    Omega_each_DS(:,:,i) = mean(matObj.Omega_est,3);
end

True_Omega = readmatrix(['./Data/GHS_sim_p',num2str(q),struc_label,'_sigmainv.csv']);

omega_elements = True_Omega(tril(true(size(True_Omega)),-1))';

tp_tn_fp_fn_matrix = zeros(length(file_index), 4);

for i =1 :length(file_index)

    temp_array = Omega_each_DS(:,:,i);
    omega_elements_current = temp_array(tril(true(size(temp_array)),-1))';
    % true Positive
   
    tp_tn_fp_fn_matrix(i,1) = sum(omega_elements~=0 & omega_elements_current~=0);
    % true negative

    tp_tn_fp_fn_matrix(i,2) = sum(omega_elements==0 & omega_elements_current==0);
    % false positive
    
    tp_tn_fp_fn_matrix(i,3) = sum(omega_elements==0 & omega_elements_current~=0);
    % false negative

    tp_tn_fp_fn_matrix(i,4) = sum(omega_elements~=0 & omega_elements_current==0);
end
MCC_matrix = zeros(1,length(file_index));

for i = 1:length(file_index)
    TP = tp_tn_fp_fn_matrix(i,1);
    TN = tp_tn_fp_fn_matrix(i,2);
    FP = tp_tn_fp_fn_matrix(i,3);
    FN = tp_tn_fp_fn_matrix(i,4);

    MCC_matrix(i) = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
end

tpr_fpr_matrix = zeros(length(file_index),2);

for i = 1: length(file_index)

    tpr_fpr_matrix(i,1) = tp_tn_fp_fn_matrix(i,1)/sum(omega_elements~=0);
    tpr_fpr_matrix(i,2) = 1- (tp_tn_fp_fn_matrix (i,2)/sum(omega_elements==0));
end

diff_Frobenious_norm = zeros(1, length(file_index));

for i = 1:length(file_index)

    diff_Frobenious_norm(i) = norm(Omega_each_DS(:,:,i) - True_Omega, 'fro');
end

stein_loss = zeros(1, length(file_index));

for i = 1:length(file_index)

    stein_loss(i) = -log(det(True_Omega\Omega_each_DS(:,:,i)))+...
        trace(True_Omega\Omega_each_DS(:,:,i))-q;

end

time_taken = zeros(1, length(file_index));

for i = 1:length(file_index)

    file_iter = file_index(i);
    FileName=['./Results/HS_LLA_',mix_type,'_mix_Workspace_of_',num2str(file_iter),'st_data_set_with_',num2str(n_EMs),...
        '_start_points_',num2str(n),'_',num2str(q),struc_label,'.mat'];

    matObj = matfile(FileName);

    time_taken(i)= mean(matObj.each_time_taken);

end


fprintf("Results for mix_type = %s\n", mix_type);
fprintf("Mean Stein Loss is,  %f, std of stein loss is  %f \n" ,round(mean(stein_loss),3), ...
    round(std(stein_loss ),3));

fprintf("Mean Fob norm is,  %f, std of From Norm is  %f \n" ,round(mean(diff_Frobenious_norm),3), ...
    round(std(diff_Frobenious_norm),3));

fprintf("Mean TPR, is,  %f, std of TPR is %f \n " ,round(mean(tpr_fpr_matrix(:,1)),3), ...
    round(std(tpr_fpr_matrix(:,1) ),3));

fprintf("Mean FPR is,  %f, std of FPR is  %f\n " ,round(mean(tpr_fpr_matrix(:,2),1),3), ...
    round(std(tpr_fpr_matrix(:,2) ),3));

fprintf("Mean MCC is,  %f, std of MCC  %f \n" ,round(mean(MCC_matrix),3), ...
    round(std(MCC_matrix),3));

fprintf("Mean Time taken is,  %f, std of Time taken is  %f \n" ,round(mean(time_taken),3), ...
    round(std(time_taken ),3));

fprintf("------------------------\n");

