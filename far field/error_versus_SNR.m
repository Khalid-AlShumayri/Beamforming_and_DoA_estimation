close all
clear all 
clc 
addpath("functions\")
%%
% Parameters
DoA = [-10 40 69];
M = 16;      % sensors
K = 3;      % number of signals (sources)
N = 100;     % number of observations (snapshots)
d = 0.5;    % Distance between elements in wavelengths
num_iterations = 100;  % Number of Monte Carlo iterations
stepsize = 1;
Res = (2/M) * 180 / pi;
L = floor(M / 2);  % Length of subarrays (spatial smoothing)


% Far-field assumption 
angles = -90:stepsize:90;  % for grid search 
a = exp(-1i * 2 * pi * d * (0:M-1)' * sin([angles(:).'] * pi / 180));
a_sps = exp(-1i * 2 * pi * d * (0:L-1)' * sin([angles(:).'] * pi / 180));
A = generate_steering_matrix(M, d, DoA);
D = generate_deriv_steering_matrix(M,d,DoA);

snr = -2:1:10;           % dB
lensnr = length(snr); 
sig_pr = 1.* ones(1, K);  % signals' power

% Define the rmse at each snr value
rmse_nor = zeros(lensnr,K);
rmse_ch  = zeros(lensnr,K);
rmse_sps = zeros(lensnr,K);

mse_nor = zeros(lensnr,K);
mse_ch  = zeros(lensnr,K); 
mse_sps = zeros(lensnr,K);

c_mu = zeros(lensnr,K);

for ir = 1:lensnr

    Pn = 10^(-snr(ir)/10);   % Noise power

    % Initialize results
    rmse_nor_mcr = zeros(num_iterations,K); % for monte carlo 
    rmse_ch_mcr  = zeros(num_iterations,K);
    rmse_sps_mcr = zeros(num_iterations,K); 

    mse_nor_mcr = zeros(num_iterations,K); % for monte carlo 
    mse_ch_mcr  = zeros(num_iterations,K);
    mse_sps_mcr = zeros(num_iterations,K); 

    C_mu_iter  = zeros(num_iterations,K);

    % Monte Carlo simulation
    tic
    for iter = 1:num_iterations
        % Generate steering matrix and signals
        
        S = diag(sig_pr./2) * (randn(K, N) + 1j * randn(K, N));
        P = S*S'./N;
        Noise = sqrt(Pn/2) * (randn(M, N) + 1j * randn(M, N));
        X = A * S + Noise;
        R = X*X'./N;
        % Generate correlated sources
        P_ch =  reshape([1 .3 .1 .3  1 .1 .1 .1 1],K,K);
    
        Lr = chol(P_ch, "lower");
        S_ch = Lr * S;
        X_ch = A * S_ch + Noise;
        R_ch = X_ch * X_ch' ./ N;
    
        % Compute spatially smoothed covariance matrix
        R_sps = zeros(L, L);
        R_ch_sps = zeros(L, L);
        for k = 1:(M - L + 1)
            x_sub = X(k:k + L - 1, :);
            x_ch = X_ch(k:k+L-1,:);
            R_sps = R_sps + (x_sub * x_sub') ./ (N + M - L + 1);
            R_ch_sps = R_ch_sps + (x_ch * x_ch') ./ (N + M - L + 1);
        end
    
        % MUSIC Algorithm for uncorrelated signals
        [Q, eigvals] = eig(R);
        [eigvals, I] = sort(diag(eigvals), 1, 'descend');
        eigvals = diag(eigvals);
        Q = Q(:, I);
        Qs = Q(:, 1:K);
        Qn = Q(:, K+1:end);
        % compute the MUSIC error variance
        H = D'*( eye(M,M)- A*inv(A'*A)*A' )*D;
        H_diag = H.*eye(K,K);
        C_mu_iter = Pn*inv(H_diag)*real(H.*(inv(P) + Pn*inv(P)*(A'*A)*inv(P)))*inv(H_diag)/(2*N);
        C_mu_iter(iter,:)= diag(C_mu_iter)';
        % ---- --- CRB -- ---- %
        
        C_crb = inv(real(H.*transpose(P)))/(2*N);
        
        % ---- ---- ---- ---- %

        % MUSIC Algorithm for correlated signals
        [Q, eigvals] = eig(R_ch);
        [eigvals, I] = sort(diag(eigvals), 1, 'descend');
        eigvals = diag(eigvals);
        Q = Q(:, I);
        Qn_ch = Q(:, K+1:end);
    
        % MUSIC Algorithm for correlated signals
        [Q, eigvals] = eig(R_ch_sps);
        [eigvals, I] = sort(diag(eigvals), 1, 'descend');
        eigvals = diag(eigvals);
        Q = Q(:, I);
        Qn_ch_sps = Q(:, K+1:L);
    
        % 1D search
        srch_music = zeros(1, length(angles));
        srch_ch    = zeros(1, length(angles));
        srch_sps= zeros(1, length(angles));
    
        for i = 1:length(angles)
            srch_music(i)  = abs(1 / (a(:, i)' * Qn * Qn' * a(:, i)));
            srch_ch(i)     = abs(1 / (a(:, i)' * Qn_ch* Qn_ch' * a(:, i)));
            srch_sps(i) = abs(1 / (a_sps(:, i)' * Qn_ch_sps * Qn_ch_sps' * a_sps(:, i)));
        end
    
        %%% ---- Peak detection ---- %%%
        [peaks, locs1] = findpeaks(srch_music);
        threshold = mean(srch_music) + 0.5*std(srch_music);
        significant_peaks = angles(locs1(peaks > threshold));
        detected_DoAs = sort(significant_peaks);
        DoA_music = zeros(1, K);
        DoA_music(1, 1:length(detected_DoAs)) = detected_DoAs;
    
        [peaks, locs2] = findpeaks(srch_ch);
        threshold = mean(srch_ch);
        significant_peaks = angles(locs2(peaks > threshold));
        detected_DoAs = sort(significant_peaks);
        DoA_ch = zeros(1, K); 
        DoA_ch(1, 1:length(detected_DoAs)) = detected_DoAs;
    
        [peaks, locs3] = findpeaks(srch_sps);
        threshold = mean(srch_sps) ;
        significant_peaks = angles(locs3(peaks > threshold));
        detected_DoAs = sort(significant_peaks);
        DoA_ch_sps = zeros(1, K);
        DoA_ch_sps(1, 1:length(detected_DoAs)) = detected_DoAs;
    
        %%% ------------------------ %%%
    
        % Calculate RMSE for this iteration
        rmse_nor_mcr(iter,:) = calc_rmse(DoA,DoA_music);
        rmse_ch_mcr(iter,:)  = calc_rmse(DoA,DoA_ch);
        rmse_sps_mcr(iter,:) = calc_rmse(DoA,DoA_ch_sps);

        mse_nor_mcr(iter,:) = calc_mse(DoA,DoA_music,d); 
        mse_ch_mcr(iter,:)  = calc_mse(DoA,DoA_ch,d);
        mse_sps_mcr(iter,:) = calc_mse(DoA,DoA_ch_sps,d);
    end

    % ----
    c_mu(ir,:) = mean(C_mu_iter);
    % Calculate and display mean RMSE
    rmse_nor(ir,:)    = mean(rmse_nor_mcr);
    rmse_ch(ir,:)     = mean(rmse_ch_mcr);
    rmse_sps(ir,:) = mean(rmse_sps_mcr);

    mse_nor(ir,:) = mean(mse_nor_mcr);
    mse_ch(ir,:) = mean(mse_ch_mcr);
    mse_sps(ir,:) = mean(mse_sps_mcr);
    
end
toc



% Plot example of the last iteration
figure;
hold on; grid on;
h1 = plot(snr, rmse_nor(:,3), 'Color', [0 .2 .9]); 
h2 = plot(snr, rmse_ch(:,3)); 
h3 = plot(snr,rmse_sps(:,3));

legend([h1 h2 h3],{'Uncorrelated', 'Correlated','Correlated + SPS'},Location="best");
title(['MUSIC algorithm : degrees']);
xlabel('SNR [dB]');
ylabel('RMSE [degrees]')

figure;
hold on; grid on;
g1 = plot(snr,mse_nor(:,3), 'Color', [0 .2 .9]); 
g2 = plot(snr,mse_ch(:,3)); 
g3 = plot(snr,mse_sps(:,3));
g4 = plot(snr,c_mu(:,3));

legend([g1 g2 g3 g4],{'Uncorrelated', 'Correlated','Correlated + SPS','MUSIC covariance'},Location="best");
title(['MUSIC algorithm : SNR vs. MSE']);
xlabel('SNR [dB]');
ylabel('MSE')

figure
plot(snr,c_mu);