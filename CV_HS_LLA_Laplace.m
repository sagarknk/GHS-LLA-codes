function CV_HS_LLA_Laplace(data_idx,n, q, prec_struc)
    
    rng(123456789);
    %%%%%%%%%%%
    matObj = matfile("HS_LLA_LAPLACE_mix.mat");
    dawson_vals = matObj.dawson_vals;
    U_grid_linear = matObj.U_grid_linear;
    %%%%%%%%%%%
    dawson_vals = dawson_vals(1:3:length(dawson_vals));
    U_grid_linear = U_grid_linear(1:3:length(U_grid_linear));
    step_size = 0.03;
    %%%%%%%%%%
    if prec_struc == 1
        struc_label = 'hubs';
    else
        struc_label = 'random';
    end
    %%%%%%%%%%
    X_mat = readmatrix(['./Data/GHS_sim_p',num2str(q),struc_label,num2str(n),'_data',num2str(data_idx),'.csv']);
    %%%%%%%%%%
    start_point = eye(q);

    for row = 2:q
        d = 0;
        while d~=q
            row_seq = row:1:q;
            col_seq = 1:1:(q-row+1);

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
    Omega_init_orig = start_point;
    fprintf("Finished start point generation\n");
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    tau_seq = exp(-3:1e-1:0);
    %tau_seq = [0.1,0.2,0.3];
    
    time_taken_taus = zeros(1, length(tau_seq));
    all_pred_MSEs = zeros(1, length(tau_seq));
    n_folds = 5;
    CV_folds = repmat(1:n_folds, 1,n/n_folds);

    for tau_idx = 1:length(tau_seq)

        pred_MSEs_this_tau = zeros(1, n_folds);
        tau_val = tau_seq(1, tau_idx);
        fprintf("Tau_idx = %d and %f tau is being used\n", tau_idx , tau_val);
        
        Omega_est = zeros(q,q,n_folds);
        total_iterations = zeros (1, n_folds);
        each_time_taken = zeros(1, n_folds);

        for CV_idx = 1:n_folds
            fprintf("%d CV-fold is being processed\n", CV_idx);
            X_mat_train = X_mat((CV_folds ~= CV_idx),:);
            X_mat_test = X_mat((CV_folds == CV_idx),:);

            S = X_mat_train' * X_mat_train;
            
            
            ind_all = zeros(q-1,q);
            for i = 1:q
                if i==1
                    ind = (2:q)';
                elseif i==q
                    ind = (1:q-1)';
                else
                    ind = [1:i-1,i+1:q]';
                end
                ind_all(:,i) = ind;
            end

            ind_all_2 = zeros(q-2,q-1);
            for i = 1:(q-1)
                if i==1
                    ind = (2:(q-1))';
                elseif i==(q-1)
                    ind = (1:q-2)';
                else
                    ind = [1:i-1,i+1:(q-1)]';
                end
                ind_all_2(:,i) = ind;
            end
            
            Omega_init = Omega_init_orig;
            Omega_current  = Omega_init_orig;
            Omega_next = eye(q);
            norm_diff = norm(Omega_current - Omega_next, 'fro');
            iter = 1;

            tic;
            while norm_diff > 1e-3

                %fprintf("%d\n",iter);
                %fprintf("%f\n", norm_diff)

                for i = 1:q
                    ind = ind_all(:,i);

                    Omega_11 = Omega_init(ind,ind);
                    s_12 = S(ind,i);
                    s_22 = S(i,i);
                    Omega_12 = Omega_init(ind,i);

                    G_B_num = zeros(1,q-1);
                    G_B_denom = zeros(1,q-1);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%% Penalty Computation %%%%%%%%%%
                    for k = 1:(q-1)

                        if(Omega_12(k,1)~=0)
                            G_B_num(1,k) = (1/tau_val)*sum(U_grid_linear...
                                .*exp(-abs(Omega_12(k,1)/tau_val).*U_grid_linear)...
                                .*dawson_vals)*step_size;

                            G_B_denom(1,k) = sum(exp(-abs(Omega_12(k,1)/tau_val)...
                                .*U_grid_linear)...
                                .*dawson_vals)*step_size;
                        else
                            G_B_num(1,k) = 1; %%% or any +ve number works
                            G_B_denom(1,k) = 0; %%% to make sure peanlty os Infinity
                        end
                    end
                    %sum(G_B_denom ==0)
                    %%%%%%%%%%%%%%%% Optimized gamma %%%%%%%%%%
                    gamma_hat = n/s_22;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    inv_Omega_11 = inv(Omega_11);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for j = 1:(q-1)

                        ind_2 = ind_all_2(:,j);
                        %C_11 = inv_Omega_11(ind_2, ind_2);
                        C_12 = inv_Omega_11(ind_2,j);
                        %C_21 = inv_Omega_11(j, ind_2);
                        C_22 = inv_Omega_11(j,j);

                        A = Omega_12(ind_2,1);
                        B = Omega_12(j,1);

                        %S_12_A = s_12(ind_2,1)';
                        s_12_B = s_12(j,1)';

                        beta_hat_num = s_12_B + s_22*(A'*C_12);
                        beta_hat_denom = s_22*C_22;

                        G_B = G_B_num(1,j)/G_B_denom(1,j);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        penalty_gamma = 2*G_B/(beta_hat_denom);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        beta_hat = beta_hat_num/beta_hat_denom;

                        if beta_hat< -1*penalty_gamma
                            Omega_12(j,1) = -penalty_gamma -beta_hat;
                        elseif beta_hat> penalty_gamma
                            Omega_12(j,1) = penalty_gamma -beta_hat;
                        else
                            Omega_12(j,1) = 0;
                        end

                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% update Omega
                    Omega_init(i,ind) = Omega_12';
                    Omega_init(ind,i) = Omega_12;
                    Omega_init(i,i) = gamma_hat + Omega_12'*inv_Omega_11* Omega_12;


                end

                Omega_next = Omega_init;
                norm_diff = norm(Omega_next - Omega_current, 'fro');
                iter = iter +1 ;
                Omega_current = Omega_init;

            end
            each_time_taken(1,CV_idx) = toc;
            total_iterations(1, CV_idx) = iter-1;
            Omega_est(:,:,CV_idx) = Omega_init;

            Coeff_matrix = Omega_init;

            for diag_idx = 1:q
                Coeff_matrix(diag_idx, diag_idx)=0;
                Coeff_matrix(:, diag_idx) = Coeff_matrix(:,diag_idx)/Omega_init(diag_idx, diag_idx);
                Coeff_matrix(:, diag_idx) = -1.*Coeff_matrix(:, diag_idx);
            end
            
            Est_matrix = X_mat_test * Coeff_matrix;
            Est_norm = norm(Est_matrix - X_mat_test, "fro");

            pred_MSEs_this_tau(1, CV_idx) = Est_norm;

        end
        
        
        %pred_MSEs_this_tau %#ok<NOPRT> 
        %mean(pred_MSEs_this_tau)
        all_pred_MSEs(1, tau_idx) = mean(pred_MSEs_this_tau);
        fprintf("%f is the avg prediction norm\n", mean(pred_MSEs_this_tau));
        time_taken_taus(1, tau_idx) = sum(each_time_taken);
    end
    
    desired_tau =  tau_seq(1, all_pred_MSEs == min(all_pred_MSEs));

    FileName = ['./Results/Tau_required_data_set_',num2str(data_idx),'_laplace_mix_',num2str(n),'_',num2str(q),struc_label,'.csv'];

    writematrix(desired_tau, FileName);
    
    fprintf("Total time taken is: %f\n", sum(time_taken_taus));
end 