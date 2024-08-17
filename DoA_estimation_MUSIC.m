close all
clear all 
clc 

% Parameters
M = 16;      % sensors
K = 2;      % number of signals (sources)
N = 50;     % number of observations (snapshots)
d = 0.5;    % Distance between elements in wavelengths
snr = 10;   %dB
snr_linear = 10^(snr/10);
sig_pr = 1 .* ones(1, K);    % signals' power
max_sig_pr = max(sig_pr);
Pn = max_sig_pr*10^(-snr/10);    % Noise power
DoA = [40 69];
num_iterations = 100;  % Number of Monte Carlo iterations
stepsize = .5;
angles = -90:stepsize:90;  % for grid search 
Res = (2/M) * 180 / pi;

L = floor(M / 2);  % Length of subarrays


% Far-field assumption 
a_sps = exp(-1i * 2 * pi * d * (0:L-1)' * sin([angles(:).'] * pi / 180));
a_coh = exp(-1i * 2 * pi * d * (0:M-1)' * sin([angles(:).'] * pi / 180));
A = generate_steering_matrix(M, d, DoA);

% Initialize results
rmse_values         = zeros(1, num_iterations);
rmse_values_ch      = zeros(1, num_iterations);
rmse_values_ch_sps  = zeros(1, num_iterations);  

% Monte Carlo simulation
for iter = 1:num_iterations
    % Generate steering matrix and signals
    
    S = diag(sqrt(sig_pr ./ 2)) * (randn(K, N) + 1j * randn(K, N));
    Noise = sqrt(Pn / 2) * (randn(M, N) + 1j * randn(M, N));
    X = A * S + Noise;

    % Generate correlated sources
    P_ch =  reshape([1 .6  1 .6],K,K);

    Lr = chol(P_ch, "lower");
    S_ch = Lr * S;
    X_ch = A * S_ch + Noise;
    R_ch = X_ch * X_ch' ./ N;

    % Compute spatially smoothed covariance matrix
    R_sps = zeros(L, L);
    R_ch_sps = zeros(L, L);;
    for k = 1:(M - L + 1)
        x_sub = X(k:k + L - 1, :);
        x_ch = X_ch(k:k+L-1,:);
        R_sps = R_sps + (x_sub * x_sub') ./ (N + M - L + 1);
        R_ch_sps = R_ch_sps + (x_ch * x_ch') ./ (N + M - L + 1);
    end

    % MUSIC Algorithm for uncorrelated signals
    [Q, D] = eig(R_sps);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qs = Q(:, 1:K);
    Qn = Q(:, K+1:L);

    % MUSIC Algorithm for correlated signals
    [Q, D] = eig(R_ch);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qn_ch = Q(:, K+1:M);

    % MUSIC Algorithm for correlated signals
    [Q, D] = eig(R_ch_sps);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qn_ch_sps = Q(:, K+1:L);


    % Grid search
    srch_sps    = zeros(1, length(angles));
    srch_ch     = zeros(1, length(angles));
    srch_ch_sps = zeros(1, length(angles));

    for i = 1:length(angles)
        srch_sps(i)    = abs(1 / (a_sps(:, i)' * Qn * Qn' * a_sps(:, i)));
        srch_ch_sps(i) = abs(1 / (a_sps(:, i)' * Qn_ch_sps * Qn_ch_sps' * a_sps(:, i)));
        srch_ch(i)     = abs(1 / (a_coh(:, i)' * Qn_ch* Qn_ch' * a_coh(:, i)));
    end

    %%% ---- Peak detection ---- %%%
    [peaks, locs1] = findpeaks(srch_sps);
    threshold = mean(srch_sps) + 0.5*std(srch_sps);
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

    [peaks, locs3] = findpeaks(srch_ch_sps);
    threshold = mean(srch_ch_sps) ;
    significant_peaks = angles(locs3(peaks > threshold));
    detected_DoAs = sort(significant_peaks);
    DoA_ch_sps = zeros(1, K);
    DoA_ch_sps(1, 1:length(detected_DoAs)) = detected_DoAs;

    %%% ------------------------ %%%

    % Calculate RMSE for this iteration
    rmse_values(iter)        = calc_rmse(DoA,DoA_music);
    rmse_values_ch(iter)     = calc_rmse(DoA,DoA_ch);
    rmse_values_ch_sps(iter) = calc_rmse(DoA,DoA_ch_sps);
end

% Calculate and display mean RMSE
rmse_sps    = mean(rmse_values);
rmse_ch     = mean(rmse_values_ch);
rmse_ch_sps = mean(rmse_values_ch_sps);

disp(['Mean RMSE (uncorrelated signals + SPS) over ',num2str(num_iterations), ' iterations: ', num2str(rmse_sps)]);
disp(['Mean RMSE (correlated signals) over ',num2str(num_iterations), ' iterations: ', num2str(rmse_ch)]);
disp(['Mean RMSE (correlated signals + SPS) over ',num2str(num_iterations), ' iterations: ', num2str(rmse_ch_sps)]);
disp(['Step-size = ',num2str(stepsize),' degrees']);

% Plot example of the last iteration
figure;
hold on; grid on;
h1 = plot(angles, srch_sps, 'Color', [0 .2 .9]); 
h2 = plot(angles, srch_ch); 
h3 = plot(angles,srch_ch_sps);
h4 = xline(DoA, '-.', 'LineWidth', .9,Color='black');
% h4 = yline(threshold, '--', 'LineWidth', .9);
legend([h1 h2 h3 h4(1)],{'Uncorrelated', 'Correlated','Correlated + SPS','Actual DoAs'},Location="best");
title(['MUSIC algorithm : SNR = ', num2str(10 * log10(min(sig_pr) / Pn)), 'dB : Resolution = ', num2str(Res), ' degrees']);
xlabel('Angle [degrees]');
ylabel('Spatial Spectrum')


%% Validation

