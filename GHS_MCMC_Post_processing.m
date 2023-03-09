n = 25;
q = 10;

prec_struc = 2; %%% can be 1 for 'hubs'

if prec_struc == 1
    struc_label = 'hubs';
else
    struc_label = 'random';
end

m = 5; %%% number of replicates 

GHS_est = zeros(q,q,m);
GHS_omega_zero = zeros(m, 0.5*q*(q-1));
GHS_Steinsloss = zeros(1,m);
GHS_Fnorm = zeros(1,m);
MCC_matrix = zeros(1,m);

sen_GHS = zeros(1,m);
spe_GHS = zeros(1,m);
prec_GHS = zeros(1,m);
acc_GHS = zeros(1,m);
fpr_GHS = zeros(1,m);
GHS_time_taken  = zeros(1,m);


for file_iter = 1:m
    fprintf("%d th data set is being processed : \n", file_iter);

    FileName=['./Results/GHS_MCMC_Workspace_of_',num2str(file_iter),'st_data_set_',num2str(n),'_',num2str(q),struc_label,'.mat'];

    matObj = matfile(FileName);
    GHS_est(:,:,file_iter) = matObj.omega_save_final;

    GHS_omega_vector_save = matObj.omega_vector_save;
    GHS_omega_vector_lb = prctile(GHS_omega_vector_save',25);
    GHS_omega_vector_ub = prctile(GHS_omega_vector_save',75);
    GHS_omega_zero(file_iter,:)= GHS_omega_vector_lb<0 & GHS_omega_vector_ub>0;
    GHS_time_taken(file_iter) = matObj.time_taken;


end

%%%% Final results Proceesing %%%%%%%%%%%%%%%%%%%%%%%%

sigma_inv = readmatrix(['./Data/GHS_sim_p',num2str(q),struc_label,'_sigmainv.csv']);
sigma = inv(sigma_inv);
omega_elements = sigma_inv(tril(true(size(sigma_inv)),-1))';    % elements in the upper triangular sigma_inv

for i = 1:m
    %fprintf("%d th data set is being processed : \n", i);
    GHS_omega_est = GHS_est(:,:,i);
    GHS_sigma_est = inv(GHS_omega_est);
    
    % Stein's loss of Sigma_inv; Frobenius norm of sigma_inv-Sigma_inv
    
    GHS_Steinsloss(i) = log(det(GHS_sigma_est*sigma_inv))+trace(GHS_omega_est*sigma)-q;
    GHS_Fnorm(i) = norm(GHS_omega_est-sigma_inv,'fro');

	TP = sum((GHS_omega_zero(i,:)==0).*(omega_elements~=0));
	FP = sum((GHS_omega_zero(i,:)==0).*(omega_elements==0));
	TN = sum((GHS_omega_zero(i,:)~=0).*(omega_elements==0));
	FN = sum((GHS_omega_zero(i,:)~=0).*(omega_elements~=0));
    
    MCC_matrix(i) = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
     
	sen_GHS(i) = TP/sum(omega_elements~=0);    % true positive rate
	spe_GHS(i) = TN/sum(omega_elements==0);
	prec_GHS(i) = TP/(TP+FP);
	acc_GHS(i) = (TP+TN)/(TP+TN+FP+FN);
	fpr_GHS(i) = 1-spe_GHS(i);
end


means = [mean(GHS_Steinsloss),mean(GHS_Fnorm),mean(sen_GHS),mean(fpr_GHS)];
fprintf('GHS means: loss, Fnorm, TPR, FPR %.3f, %.3f, %.3f, %.3f\n',means)

sds = [std(GHS_Steinsloss),std(GHS_Fnorm),std(sen_GHS),std(fpr_GHS)];
fprintf('GHS sds: loss, Fnorm, TPR, FPR %.3f, %.3f, %.3f, %.3f\n',sds)

fprintf("Mean MCC is,  %f, std of MCC  %f \n" ,round(mean(MCC_matrix),3), ...
    round(std(MCC_matrix),3));

avg_time_taken = mean(GHS_time_taken);
fprintf('%.3f is the avg time taken\n', avg_time_taken);
 